using IQGenetics.MC1R.Core.Alignment;
using IQGenetics.MC1R.Core.Sanger;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace IQGenetics.MC1R.Core.MC1R;

public sealed record Mc1rRead2026Result(
    string SampleName,
    string FilePath,
    Orientation Orientation,
    int AlignmentScore,
    bool IsDirty,
    string DirtyReason,
    SiteCall C212,
    SiteCall C274,
    SiteCall C376,
    SiteCall C398,
    SiteCall C409,
    SiteCall C427,
    SiteCall C637,
    SiteCall C644,
    ReadQcSummary QcSummary
);

/// <summary>
/// 2026 MC1R per-read caller. Reuses the validated global chromatogram QC from
/// the legacy caller but calls the eight missense sites used by Ma et al. 2026.
/// </summary>
public static class Mc1r2026Caller
{
    public const int SiteC212 = 212;
    public const int SiteC274 = 274;
    public const int SiteC376 = 376;
    public const int SiteC398 = 398;
    public const int SiteC409 = 409;
    public const int SiteC427 = 427;
    public const int SiteC637 = 637;
    public const int SiteC644 = 644;

    public static Mc1rRead2026Result CallAb1(Ab1Chromatogram ab1, Mc1rReference reference)
    {
        var (trimStart, trimEnd) = FindTrim(ab1.Qualities, qmin: 15);
        if (trimEnd - trimStart < 150)
            throw new InvalidDataException("AB1 read too short after trimming.");

        var query = ab1.Bases.Substring(trimStart, trimEnd - trimStart);
        var alnF = SmithWaterman.Align(reference.CdsSequence, query);
        var alnRC = SmithWaterman.Align(reference.CdsSequence, ReverseComplement(query));
        Orientation orient = alnF.Score >= alnRC.Score ? Orientation.Forward : Orientation.ReverseComplement;
        var aln = orient == Orientation.Forward ? alnF : alnRC;

        if (aln.Score < 600)
            return NoCall(ab1, orient, aln.Score, "Low alignment score — reference mismatch or very poor sequencing.");

        var qc = Mc1rCaller.EvaluateReadQc(ab1, trimStart, query.Length, aln, orient);
        if (qc.IsDirty)
            return NoCall(ab1, orient, aln.Score, qc.Reason, qc);

        return new Mc1rRead2026Result(
            ab1.FileName, ab1.FilePath, orient, aln.Score, false, "",
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC212),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC274),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC376),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC398),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC409),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC427),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC637),
            CallSite(ab1, trimStart, query.Length, aln, orient, SiteC644),
            qc
        );
    }

    private static Mc1rRead2026Result NoCall(Ab1Chromatogram ab1, Orientation orient, int score, string reason, ReadQcSummary? qc = null)
    {
        qc ??= new ReadQcSummary(true, reason, 0, 0, 0, 0, 0, 0, 0, 0, "Mixed template (persistent)");
        SiteCall N(int p) => new(p, "NoCall", reason);
        return new Mc1rRead2026Result(
            ab1.FileName, ab1.FilePath, orient, score, true, reason,
            N(212), N(274), N(376), N(398), N(409), N(427), N(637), N(644), qc
        );
    }

    private static (int Start, int End) FindTrim(int[] qualities, int qmin)
    {
        int start = 0;
        while (start < qualities.Length && qualities[start] < qmin) start++;
        int end = qualities.Length;
        while (end > start && qualities[end - 1] < qmin) end--;
        return (start, end);
    }

    private static SiteCall CallSite(Ab1Chromatogram ab1, int trimStart, int trimmedLen, SmithWaterman.Result aln, Orientation orient, int cdnaPos)
    {
        if (!aln.RefToQueryMap.TryGetValue(cdnaPos - 1, out int qIdxTransformed))
            return new SiteCall(cdnaPos, "NoCall", "Reference position not mapped.");

        int qIdxOriginal = orient == Orientation.Forward ? qIdxTransformed : (trimmedLen - 1 - qIdxTransformed);
        int bidx = trimStart + qIdxOriginal;
        return GenotypeAt(ab1, bidx, orient, cdnaPos);
    }

    private sealed record PeakSummary(char TopBase, char SecondBase, int Top, int Second, int Sum, double SecondFraction, double SecondOverTop);

    private static SiteCall GenotypeAt(Ab1Chromatogram ab1, int bidx, Orientation orient, int cdnaPos)
    {
        const int qMin = 15;
        if (bidx < 0 || bidx >= ab1.Bases.Length)
            return new SiteCall(cdnaPos, "NoCall", "Index out of range.");

        var peak = PeakAtBase(ab1, bidx, 2);
        if (peak.Sum < 200)
            return new SiteCall(cdnaPos, "NoCall", "Low signal at site.");

        bool isHet = peak.SecondFraction >= 0.22 && peak.SecondOverTop >= 0.33;
        int q = ab1.Qualities[bidx];
        string? note = null;

        if (q < qMin)
        {
            // Do not lower the global quality threshold. A low-Q site is rescued only when
            // the entire chromatogram has already passed QC and the site shows a strong,
            // localized dual peak consistent with genuine Sanger heterozygosity.
            if (!(q >= 8 && isHet && peak.Sum >= 400 && peak.SecondOverTop >= 0.45))
                return new SiteCall(cdnaPos, "NoCall", $"Low Q at site (Q={q}).");
            note = $"Low-Q strong dual peak accepted after global QC (Q={q}).";
        }

        char a1 = peak.TopBase;
        char a2 = isHet ? peak.SecondBase : peak.TopBase;
        if (orient == Orientation.ReverseComplement)
        {
            a1 = Complement(a1);
            a2 = Complement(a2);
        }

        var alleles = new[] { a1, a2 };
        Array.Sort(alleles);
        return new SiteCall(cdnaPos, $"{alleles[0]}/{alleles[1]}", note);
    }

    private static PeakSummary PeakAtBase(Ab1Chromatogram ab1, int baseIndex, int window)
    {
        int x = ab1.PeakLocations[baseIndex];
        var values = new List<(char Base, int Val)>
        {
            ('A', MaxInWindow(ab1, 'A', x, window)),
            ('C', MaxInWindow(ab1, 'C', x, window)),
            ('G', MaxInWindow(ab1, 'G', x, window)),
            ('T', MaxInWindow(ab1, 'T', x, window))
        };
        values.Sort((u, v) => v.Val.CompareTo(u.Val));
        int sum = values.Sum(v => v.Val);
        var top = values[0];
        var second = values[1];
        double secFrac = sum > 0 ? (double)second.Val / sum : 0;
        double secOverTop = top.Val > 0 ? (double)second.Val / top.Val : 0;
        return new PeakSummary(top.Base, second.Base, top.Val, second.Val, sum, secFrac, secOverTop);
    }

    private static int MaxInWindow(Ab1Chromatogram ab1, char baseChar, int x, int w)
    {
        if (!ab1.Traces.TryGetValue(baseChar, out var trace)) return 0;
        int lo = Math.Max(0, x - w);
        int hi = Math.Min(trace.Length - 1, x + w);
        short max = 0;
        for (int i = lo; i <= hi; i++) if (trace[i] > max) max = trace[i];
        return max;
    }

    private static string ReverseComplement(string s)
    {
        var arr = s.ToCharArray();
        Array.Reverse(arr);
        for (int i = 0; i < arr.Length; i++) arr[i] = Complement(arr[i]);
        return new string(arr);
    }

    private static char Complement(char b) => b switch
    {
        'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N'
    };
}
