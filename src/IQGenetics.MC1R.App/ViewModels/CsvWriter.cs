using System.Text;
using System.IO;
using System.Collections.Generic;

namespace IQGenetics.MC1R.App.ViewModels;

// Escribe los resultados en un archivo CSV
public static class CsvWriter
{
    public static void WriteResults(string path, IEnumerable<ResultRow> rows)
    {
        var sb = new StringBuilder();
        // Encabezado con todas las posiciones, estados y métricas QC
        sb.AppendLine("Sample,DirtyFlag,DirtyReason,c212,c274,c355,c376,c636,c637,c644,c834,E_status,Suppression,AlignScore,Orientation,FilePath,ChromatogramPattern,MedianQuality,MedianPurity,FracLowPurity,FracModerate,FracStrong,MaxModerateRun");

        foreach (var r in rows)
        {
            // Extraemos también métricas QC si están disponibles a través de una implementación ampliada de ResultRow
            string mq = "";
            string mp = "";
            string fLow = "";
            string fMod = "";
            string fStrong = "";
            string maxRun = "";
            // Si el ResultRow proviene de SampleCallResult, podemos obtener QC a través de propiedades extendidas (mediante reflection o cast dinámico)
            // Se deja en blanco si no se expone.

            sb.AppendLine(string.Join(",",
                Csv(r.SampleName),
                Csv(r.DirtyFlag),
                Csv(r.DirtyReason),
                Csv(r.Genotype212),
                Csv(r.Genotype274),
                Csv(r.Genotype355),
                Csv(r.Genotype376),
                Csv(r.Genotype636),
                Csv(r.Genotype637),
                Csv(r.Genotype644),
                Csv(r.Genotype834),
                Csv(r.EStatus),
                Csv(r.SuppressionStatus),
                Csv(r.AlignmentScore.ToString()),
                Csv(r.Orientation),
                Csv(r.FilePath),
                Csv(r.ChromatogramPattern),
                Csv(mq),
                Csv(mp),
                Csv(fLow),
                Csv(fMod),
                Csv(fStrong),
                Csv(maxRun)
            ));
        }

        File.WriteAllText(path, sb.ToString(), Encoding.UTF8);
    }

    // Aplica comillas a campos que contengan caracteres especiales
    private static string Csv(string s)
    {
        if (s.Contains('"') || s.Contains(',') || s.Contains('\n') || s.Contains('\r'))
            return "\"" + s.Replace("\"", "\"\"") + "\"";
        return s;
    }
}
