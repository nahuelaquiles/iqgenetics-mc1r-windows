using System.Text;
using System.IO;
using System.Collections.Generic;

namespace IQGenetics.MC1R.App.ViewModels
{
    public static class CsvWriter
    {
        public static void WriteResults(string path, IEnumerable<ResultRow> rows)
        {
            var sb = new StringBuilder();
            sb.AppendLine("Sample,TotalReads,UsableReads,QC,QCNotes,c212,c274,c376,c398,c409,c427,c637,c644,CallStatus,BreedingCategory,CompatibleDiplotypes,CompatiblePhenotypes,Interpretation");
            foreach (var r in rows)
            {
                sb.AppendLine(string.Join(",",
                    Csv(r.SampleName), Csv(r.TotalReads.ToString()), Csv(r.UsableReads.ToString()), Csv(r.QcStatus), Csv(r.QcNotes),
                    Csv(r.Genotype212), Csv(r.Genotype274), Csv(r.Genotype376), Csv(r.Genotype398), Csv(r.Genotype409), Csv(r.Genotype427), Csv(r.Genotype637), Csv(r.Genotype644),
                    Csv(r.CallStatus), Csv(r.BreedingCategory), Csv(r.CompatibleDiplotypes), Csv(r.CompatiblePhenotypes), Csv(r.Interpretation)));
            }
            File.WriteAllText(path, sb.ToString(), new UTF8Encoding(true));
        }

        private static string Csv(string s)
        {
            if (s.Contains('"') || s.Contains(',') || s.Contains('\n') || s.Contains('\r') || s.Contains(';'))
                return "\"" + s.Replace("\"", "\"\"") + "\"";
            return s;
        }
    }
}
