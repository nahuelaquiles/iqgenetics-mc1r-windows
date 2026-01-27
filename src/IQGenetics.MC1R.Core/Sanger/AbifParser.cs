using System.Buffers.Binary;
using System.Text;

namespace IQGenetics.MC1R.Core.Sanger;

/// <summary>
/// Minimal ABIF (AB1) parser for Sanger chromatograms.
/// Extracts: PBAS1 (basecalls), PCON1 (quality bytes), PLOC1 (peak locations), FWO_1 (channel order), DATA9-12 (trace channels).
/// </summary>
public static class AbifParser
{
    public static Ab1Chromatogram ReadAb1(string path)
    {
        var bytes = File.ReadAllBytes(path);
        if (bytes.Length < 64) throw new InvalidDataException("File too small.");
        if (Encoding.ASCII.GetString(bytes, 0, 4) != "ABIF") throw new InvalidDataException("Not an AB1/ABIF file.");

        // Root directory entry begins at offset 6 (ABIF spec)
        var root = ReadEntry(bytes, 6);
        var dirOffset = (int)root.DataOffset;
        var count = (int)root.NumElements;

        if (dirOffset <= 0 || dirOffset + count * 28 > bytes.Length)
            throw new InvalidDataException("Invalid ABIF directory.");

        var entries = new List<AbifDirectoryEntry>(count);
        for (int i = 0; i < count; i++)
        {
            entries.Add(ReadEntry(bytes, dirOffset + i * 28));
        }

        byte[] pbas = GetRaw(entries, bytes, "PBAS", 1);
        byte[] pcon = GetRaw(entries, bytes, "PCON", 1);
        short[] ploc = ReadShortArray(GetRaw(entries, bytes, "PLOC", 1));

        string fwo = ReadAscii(GetRaw(entries, bytes, "FWO_", 1)).Trim();
        if (fwo.Length != 4) fwo = "GATC"; // safe default

        // Channels are usually DATA9..DATA12 in order given by FWO_1
        var data9 = ReadShortArray(GetRaw(entries, bytes, "DATA", 9));
        var data10 = ReadShortArray(GetRaw(entries, bytes, "DATA", 10));
        var data11 = ReadShortArray(GetRaw(entries, bytes, "DATA", 11));
        var data12 = ReadShortArray(GetRaw(entries, bytes, "DATA", 12));

        if (pbas.Length != pcon.Length || pbas.Length != ploc.Length)
            throw new InvalidDataException("Inconsistent PBAS/PCON/PLOC lengths.");

        var traces = new Dictionary<char, short[]>();
        // Map FWO letters to the corresponding DATA arrays
        var datas = new[] { data9, data10, data11, data12 };
        for (int i = 0; i < 4; i++)
        {
            traces[fwo[i]] = datas[i];
        }

        string bases = ReadAscii(pbas);
        var quals = pcon.Select(b => (int)b).ToArray();

        return new Ab1Chromatogram(path, bases, quals, ploc, fwo, traces);
    }

    private static AbifDirectoryEntry ReadEntry(byte[] bytes, int offset)
    {
        string name = Encoding.ASCII.GetString(bytes, offset, 4);
        uint number = ReadUInt32BE(bytes, offset + 4);
        var etype = (AbifElementType)ReadUInt16BE(bytes, offset + 8);
        ushort esize = ReadUInt16BE(bytes, offset + 10);
        uint n = ReadUInt32BE(bytes, offset + 12);
        uint dsize = ReadUInt32BE(bytes, offset + 16);
        uint doff = ReadUInt32BE(bytes, offset + 20);
        uint handle = ReadUInt32BE(bytes, offset + 24);
        return new AbifDirectoryEntry(name, number, etype, esize, n, dsize, doff, handle);
    }

    private static byte[] GetRaw(List<AbifDirectoryEntry> entries, byte[] bytes, string name, uint number)
    {
        var e = entries.FirstOrDefault(x => x.Name == name && x.Number == number);
        if (e is null) throw new InvalidDataException($"Missing AB1 tag {name}{number}.");
        // If DataSize <= 4, DataOffset stores the data inline.
        if (e.DataSize <= 4)
        {
            // DataOffset is a 4-byte big-endian value. The first DataSize bytes represent the inline data.
            Span<byte> tmp = stackalloc byte[4];
            BinaryPrimitives.WriteUInt32BigEndian(tmp, e.DataOffset);
            return tmp[..(int)e.DataSize].ToArray();
        }

        int start = (int)e.DataOffset;
        int len = (int)e.DataSize;
        if (start < 0 || start + len > bytes.Length) throw new InvalidDataException($"Invalid offset for {name}{number}.");
        var raw = new byte[len];
        Buffer.BlockCopy(bytes, start, raw, 0, len);
        return raw;
    }

    private static string ReadAscii(byte[] raw) => Encoding.ASCII.GetString(raw);

    private static ushort ReadUInt16BE(byte[] b, int offset) => BinaryPrimitives.ReadUInt16BigEndian(b.AsSpan(offset, 2));
    private static uint ReadUInt32BE(byte[] b, int offset) => BinaryPrimitives.ReadUInt32BigEndian(b.AsSpan(offset, 4));

    private static short[] ReadShortArray(byte[] raw)
    {
        if (raw.Length % 2 != 0) throw new InvalidDataException("Short array byte length must be even.");
        var arr = new short[raw.Length / 2];
        for (int i = 0; i < arr.Length; i++)
        {
            arr[i] = BinaryPrimitives.ReadInt16BigEndian(raw.AsSpan(i * 2, 2));
        }
        return arr;
    }
}
