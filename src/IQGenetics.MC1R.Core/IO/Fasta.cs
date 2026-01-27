using System.Text;
using System.IO;


namespace IQGenetics.MC1R.Core.IO;

public sealed record FastaRecord(string Header, string Sequence);

public static class Fasta
{
    public static FastaRecord ReadSingle(string path)
    {
        var lines = File.ReadAllLines(path);
        if (lines.Length == 0) throw new InvalidDataException("Empty FASTA.");
        string header = "";
        var sb = new StringBuilder();

        foreach (var raw in lines)
        {
            var line = raw.Trim();
            if (line.Length == 0) continue;
            if (line.StartsWith(">"))
            {
                if (header.Length > 0)
                    throw new InvalidDataException("FASTA contains multiple records. Provide a single reference record.");
                header = line[1..].Trim();
            }
            else
            {
                sb.Append(line);
            }
        }

        if (sb.Length == 0) throw new InvalidDataException("No sequence found in FASTA.");
        var seq = Normalize(sb.ToString());
        return new FastaRecord(header, seq);
    }

    public static string Normalize(string seq)
    {
        var sb = new StringBuilder(seq.Length);
        foreach (var c in seq)
        {
            char u = char.ToUpperInvariant(c);
            if (u is 'A' or 'C' or 'G' or 'T' or 'N')
                sb.Append(u);
        }
        return sb.ToString();
    }
}
