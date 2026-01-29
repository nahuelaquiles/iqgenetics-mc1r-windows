using IQGenetics.MC1R.Core.Alignment;
using IQGenetics.MC1R.Core.Sanger;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace IQGenetics.MC1R.Core.MC1R;

// Representa la orientación de la lectura
public enum Orientation
{
    Forward,
    ReverseComplement
}

// Contiene la llamada a un sitio concreto (posición cDNA, genotipo y nota opcional)
public sealed record SiteCall(
    int CdnaPosition,
    string Genotype,
    string? Notes
);

// Contiene todas las llamadas de un archivo de entrada
public sealed record SampleCallResult(
    string SampleName,
    string FilePath,
    Orientation Orientation,
    int AlignmentScore,
    bool IsDirty,
    string DirtyReason,
    SiteCall C212,
    SiteCall C274,
    SiteCall C355,
    SiteCall C376,
    SiteCall C636,
    SiteCall C637,
    SiteCall C644,
    SiteCall C834,
    string EStatus,
    string SuppressionStatus
);

public static class Mc1rCaller
{
    // Sitios de interés según el SOP (c.212, c.274, c.355, c.376, c.636, c.637, c.644 y c.834)
    public const int SiteC212 = 212;
    public const int SiteC274 = 274;
    public const int SiteC355 = 355;
    public const int SiteC376 = 376;
    public const int SiteC636 = 636;
    public const int SiteC637 = 637;
    public const int SiteC644 = 644;
    public const int SiteC834 = 834;

    public static SampleCallResult CallAb1(Ab1Chromatogram ab1, Mc1rReference reference)
    {
        // Recorta los extremos de mala calidad pero mantiene la correspondencia de índices
        var (trimStart, trimEnd) = FindTrim(ab1.Qualities, qmin: 15);
        if (trimEnd - trimStart < 150)
            throw new InvalidDataException("AB1 read too short after trimming.");

        var query = ab1.Bases.Substring(trimStart, trimEnd - trimStart);
        var alnF = SmithWaterman.Align(reference.CdsSequence, query);
        var alnRC = SmithWaterman.Align(reference.CdsSequence, ReverseComplement(query));

        Orientation orient = alnF.Score >= alnRC.Score ? Orientation.Forward : Orientation.ReverseComplement;
        var aln = orient == Orientation.Forward ? alnF : alnRC;

        // Comprueba que el alineamiento sea suficientemente bueno
        if (aln.Score < 600)
        {
            return new SampleCallResult(
                ab1.FileName, ab1.FilePath, orient, aln.Score,
                true, "Low alignment score — reference mismatch or very poor sequencing.",
                new SiteCall(SiteC212, "NoCall", "Alignment failed"),
                new SiteCall(SiteC274, "NoCall", "Alignment failed"),
                new SiteCall(SiteC355, "NoCall", "Alignment failed"),
                new SiteCall(SiteC376, "NoCall", "Alignment failed"),
                new SiteCall(SiteC636, "NoCall", "Alignment failed"),
                new SiteCall(SiteC637, "NoCall", "Alignment failed"),
                new SiteCall(SiteC644, "NoCall", "Alignment failed"),
                new SiteCall(SiteC834, "NoCall", "Alignment failed"),
                "NoCall", "NoCall"
            );
        }

        // Detección de suciedad (mezcla de templados o baja pureza de picos)
        var dirty = ComputeDirtyMetrics(ab1, trimStart, query.Length, aln, orient);
        if (dirty.IsDirty)
        {
            return new SampleCallResult(
                ab1.FileName, ab1.FilePath, orient, aln.Score,
                true, dirty.Reason,
                new SiteCall(SiteC212, "NoCall", "DIRTY"),
                new SiteCall(SiteC274, "NoCall", "DIRTY"),
                new SiteCall(SiteC355, "NoCall", "DIRTY"),
                new SiteCall(SiteC376, "NoCall", "DIRTY"),
                new SiteCall(SiteC636, "NoCall", "DIRTY"),
                new SiteCall(SiteC637, "NoCall", "DIRTY"),
                new SiteCall(SiteC644, "NoCall", "DIRTY"),
                new SiteCall(SiteC834, "NoCall", "DIRTY"),
                "DIRTY", "DIRTY"
            );
        }

        // Llama a todos los sitios para obtener sus genotipos
        var c212 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC212);
        var c274 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC274);
        var c355 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC355);
        var c376 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC376);
        var c636 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC636);
        var c637 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC637);
        var c644 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC644);
        var c834 = CallSite(ab1, trimStart, query.Length, aln, orient, SiteC834);

        // Interpreta los alelos principales
        string eStatus = InterpretE(c274.Genotype);
        string suppression = InterpretSuppression(c644.Genotype);

        return new SampleCallResult(
            ab1.FileName, ab1.FilePath, orient, aln.Score,
            false, "",
            c212, c274, c355, c376, c636, c637, c644, c834,
            eStatus, suppression
        );
    }

    // Determina el segmento utilizable de la secuencia según la calidad mínima
    private static (int Start, int End) FindTrim(int[] qualities, int qmin)
    {
        int start = 0;
        while (start < qualities.Length && qualities[start] < qmin) start++;
        int end = qualities.Length;
        while (end > start && qualities[end - 1] < qmin) end--;
        return (start, end);
    }

    private sealed record DirtyMetrics(bool IsDirty, string Reason);

    // Calcula si la muestra está sucia en función de la pureza de los picos
    private static DirtyMetrics ComputeDirtyMetrics(Ab1Chromatogram ab1, int trimStart, int trimmedLen, SmithWaterman.Result aln, Orientation orient)
    {
        const int Qmin = 20;
        const int window = 2;
        const int minSum = 200;

        var purities = new List<double>(capacity: 800);
        int total = 0;
        int lowPurity = 0;

        foreach (var kv in aln.RefToQueryMap)
        {
            int refIdx = kv.Key;
            int qIdxTransformed = kv.Value;
            int qIdxOriginal = orient == Orientation.Forward
                ? qIdxTransformed
                : (trimmedLen - 1 - qIdxTransformed);

            int bidx = trimStart + qIdxOriginal;
            if (bidx < 0 || bidx >= ab1.Qualities.Length) continue;
            if (ab1.Qualities[bidx] < Qmin) continue;

            var peak = PeakAtBase(ab1, bidx, window);
            if (peak.Sum < minSum) continue;

            total++;
            purities.Add(peak.Purity);
            if (peak.Purity < 0.55) lowPurity++;
        }

        if (total < 200)
            return new DirtyMetrics(true, "Insufficient high-quality aligned region.");

        double fracLow = (double)lowPurity / total;
        double medianPurity = Median(purities);
        double medianQ = Median(ab1.Qualities.Select(q => (double)q).ToList());

        if (medianQ < 20)
            return new DirtyMetrics(true, "Low base-call quality across read.");

        if (fracLow > 0.20 || medianPurity < 0.65)
            return new DirtyMetrics(true, "DIRTY / MIXED TEMPLATE — repeat PCR/sequencing (persistent secondary peaks).");

        return new DirtyMetrics(false, "");
    }

    // Llama a un sitio concreto utilizando el mapa de alineamiento
    private static SiteCall CallSite(Ab1Chromatogram ab1, int trimStart, int trimmedLen, SmithWaterman.Result aln, Orientation orient, int cdnaPos)
    {
        int refIdx = cdnaPos - 1;
        if (!aln.RefToQueryMap.TryGetValue(refIdx, out int qIdxTransformed))
            return new SiteCall(cdnaPos, "NoCall", "Reference position not mapped.");

        int qIdxOriginal = orient == Orientation.Forward ? qIdxTransformed : (trimmedLen - 1 - qIdxTransformed);
        int bidx = trimStart + qIdxOriginal;

        var call = GenotypeAt(ab1, bidx, orient);
        return new SiteCall(cdnaPos, call.Genotype, call.Notes);
    }

    private sealed record GenotypeCall(string Genotype, string? Notes);

    // Determina el genotipo en un índice determinado valorando heterocigosis
    private static GenotypeCall GenotypeAt(Ab1Chromatogram ab1, int bidx, Orientation orient)
    {
        const int Qmin = 15;
        const int window = 2;
        const int minSum = 200;

        if (bidx < 0 || bidx >= ab1.Bases.Length) return new GenotypeCall("NoCall", "Index out of range.");
        if (ab1.Qualities[bidx] < Qmin) return new GenotypeCall("NoCall", "Low Q at site.");

        var peak = PeakAtBase(ab1, bidx, window);
        if (peak.Sum < minSum) return new GenotypeCall("NoCall", "Low signal at site.");

        bool isHet = peak.SecondFraction >= 0.22 && peak.SecondOverTop >= 0.33;

        char a1 = peak.TopBase;
        char a2 = isHet ? peak.SecondBase : peak.TopBase;

        if (orient == Orientation.ReverseComplement)
        {
            a1 = Complement(a1);
            a2 = Complement(a2);
        }

        var alleles = new[] { a1, a2 };
        Array.Sort(alleles);
        string genotype = $"{alleles[0]}/{alleles[1]}";

        return new GenotypeCall(genotype, null);
    }

    private sealed record PeakSummary(char TopBase, char SecondBase, int Top, int Second, int Sum, double Purity, double SecondFraction, double SecondOverTop);

    // Calcula las intensidades máximas en una ventana para cada base
    private static PeakSummary PeakAtBase(Ab1Chromatogram ab1, int baseIndex, int window)
    {
        int x = ab1.PeakLocations[baseIndex];

        int A = MaxInWindow(ab1, 'A', x, window);
        int C = MaxInWindow(ab1, 'C', x, window);
        int G = MaxInWindow(ab1, 'G', x, window);
        int T = MaxInWindow(ab1, 'T', x, window);

        var list = new List<(char Base, int Val)>
        {
            ('A', A), ('C', C), ('G', G), ('T', T)
        };
        list.Sort((u, v) => v.Val.CompareTo(u.Val));

        var top = list[0];
        var second = list[1];

        int sum = A + C + G + T;
        double purity = sum > 0 ? (double)top.Val / sum : 0.0;
        double secFrac = sum > 0 ? (double)second.Val / sum : 0.0;
        double secOverTop = top.Val > 0 ? (double)second.Val / top.Val : 0.0;

        return new PeakSummary(top.Base, second.Base, top.Val, second.Val, sum, purity, secFrac, secOverTop);
    }

    private static int MaxInWindow(Ab1Chromatogram ab1, char baseChar, int x, int w)
    {
        if (!ab1.Traces.TryGetValue(baseChar, out var trace)) return 0;
        int lo = Math.Max(0, x - w);
        int hi = Math.Min(trace.Length - 1, x + w);
        short max = 0;
        for (int i = lo; i <= hi; i++)
        {
            if (trace[i] > max) max = trace[i];
        }
        return max;
    }

    // Interpreta el alelo E según el genotipo en c.274
    private static string InterpretE(string genotype274)
    {
        if (genotype274 == "NoCall") return "NoCall";
        bool hasA = genotype274.Contains('A');
        if (!hasA) return "Not-E";
        if (genotype274 == "A/A") return "E/E";
        return "E/other";
    }

    // Interpreta el estado de supresión según el genotipo en c.644
    private static string InterpretSuppression(string genotype644)
    {
        if (genotype644 == "NoCall") return "NoCall";
        return genotype644.Contains('C') ? "Suppressed" : "Not suppressed";
    }

    // Calcula la cadena complementaria inversa de una secuencia
    private static string ReverseComplement(string s)
    {
        var arr = s.ToCharArray();
        Array.Reverse(arr);
        for (int i = 0; i < arr.Length; i++) arr[i] = Complement(arr[i]);
        return new string(arr);
    }

    // Devuelve el complemento de una base
    private static char Complement(char b) => b switch
    {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => 'N'
    };

    // Calcula la mediana de una lista de valores
    private static double Median(List<double> vals)
    {
        if (vals.Count == 0) return 0;
        vals.Sort();
        int mid = vals.Count / 2;
        if (vals.Count % 2 == 1) return vals[mid];
        return 0.5 * (vals[mid - 1] + vals[mid]);
    }
}
