using System.Text;

namespace IQGenetics.MC1R.App.ViewModels;

public static class CsvWriter
{
    public static void WriteResults(string path, IEnumerable<ResultRow> rows)
    {
        var sb = new StringBuilder();
        sb.AppendLine("Sample,DirtyFlag,DirtyReason,c274,E_status,c644,Suppression,AlignScore,Orientation,FilePath");

        foreach (var r in rows)
        {
            sb.AppendLine(string.Join(",",
                Csv(r.SampleName),
                Csv(r.DirtyFlag),
                Csv(r.DirtyReason),
                Csv(r.Genotype274),
                Csv(r.EStatus),
                Csv(r.Genotype644),
                Csv(r.SuppressionStatus),
                Csv(r.AlignmentScore.ToString()),
                Csv(r.Orientation),
                Csv(r.FilePath)
            ));
        }

        File.WriteAllText(path, sb.ToString(), Encoding.UTF8);
    }

    private static string Csv(string s)
    {
        if (s.Contains('"') || s.Contains(',') || s.Contains('\n') || s.Contains('\r'))
            return "\"" + s.Replace("\"", "\"\"") + "\"";
        return s;
    }
}
