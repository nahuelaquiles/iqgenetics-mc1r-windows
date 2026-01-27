using IQGenetics.MC1R.Core.IO;

namespace IQGenetics.MC1R.Core.MC1R;

public sealed record Mc1rReference(string Header, string CdsSequence, string SourcePath);

public static class Mc1rReferenceLoader
{
    /// <summary>
    /// Loads a reference file. If the sequence is longer than typical CDS, attempts to extract a 945 bp ORF (MC1R CDS length).
    /// Assumes the sequence is in the correct orientation (most NCBI exports are).
    /// </summary>
    public static Mc1rReference Load(string path)
    {
        var ext = Path.GetExtension(path).ToLowerInvariant();
        string header = "";
        string seq;

        if (ext is ".fa" or ".fasta" or ".fna" or ".fas" or ".txt")
        {
            var rec = Fasta.ReadSingle(path);
            header = rec.Header;
            seq = rec.Sequence;
        }
        else
        {
            // fallback: treat as plain text sequence
            seq = SequenceText.ReadPlain(path);
        }

        if (seq.Length < 200)
            throw new InvalidDataException("Reference sequence too short.");

        // If already looks like CDS (starts with ATG and ~945 bp), accept.
        if (seq.StartsWith("ATG") && seq.Length is >= 900 and <= 1100)
        {
            return new Mc1rReference(header, seq, path);
        }

        // Try to find a 945 bp ORF starting with ATG and ending with a stop (TAG/TAA/TGA).
        var cds = TryExtractOrf(seq, orfLength: 945);
        if (cds is null)
        {
            throw new InvalidDataException("Could not locate a 945 bp MC1R CDS ORF in the provided reference. Please provide CDS-only reference starting at ATG.");
        }

        return new Mc1rReference(header, cds, path);
    }

    private static string? TryExtractOrf(string seq, int orfLength)
    {
        string[] stops = ["TAA", "TAG", "TGA"];
        // scan all frames
        for (int frame = 0; frame < 3; frame++)
        {
            for (int i = frame; i + 3 <= seq.Length; i += 3)
            {
                if (seq.AsSpan(i, 3).SequenceEqual("ATG".AsSpan()))
                {
                    int j = i + orfLength - 3; // last codon start
                    if (j + 3 > seq.Length) continue;
                    var stop = seq.Substring(j, 3);
                    if (stops.Contains(stop))
                    {
                        return seq.Substring(i, orfLength);
                    }
                }
            }
        }
        return null;
    }
}
