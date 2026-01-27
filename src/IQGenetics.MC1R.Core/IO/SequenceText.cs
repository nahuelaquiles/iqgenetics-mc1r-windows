using System.IO;
namespace IQGenetics.MC1R.Core.IO;

public static class SequenceText
{
    /// <summary>Reads a plain sequence export (SEQ/TXT) and keeps only A/C/G/T/N.</summary>
    public static string ReadPlain(string path)
    {
        var text = File.ReadAllText(path);
        return Fasta.Normalize(text);
    }
}
