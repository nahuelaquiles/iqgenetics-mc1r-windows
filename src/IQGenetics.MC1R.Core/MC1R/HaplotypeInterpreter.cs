using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace IQGenetics.MC1R.Core.MC1R;

public sealed record Mc1rHaplotype(
    string Id,
    string Phenotype,
    char C212,
    char C274,
    char C376,
    char C398,
    char C409,
    char C427,
    char C637,
    char C644
);

public sealed record ConsensusSiteCall(
    int CdnaPosition,
    string Genotype,
    int SupportingReads,
    bool IsDiscordant,
    string Note
);

public sealed record SampleConsensusResult(
    string SampleId,
    int TotalReads,
    int UsableReads,
    string QcStatus,
    string QcNotes,
    ConsensusSiteCall C212,
    ConsensusSiteCall C274,
    ConsensusSiteCall C376,
    ConsensusSiteCall C398,
    ConsensusSiteCall C409,
    ConsensusSiteCall C427,
    ConsensusSiteCall C637,
    ConsensusSiteCall C644,
    string CallStatus,
    string BreedingCategory,
    string CompatibleDiplotypes,
    string CompatiblePhenotypes,
    string Interpretation
);

public static class Mc1rHaplotypeInterpreter
{
    public static IReadOnlyList<Mc1rHaplotype> CommonHaplotypes { get; } = new[]
    {
        new Mc1rHaplotype("N1",   "Wild-type",               'T','G','G','T','G','A','C','A'),
        new Mc1rHaplotype("N2",   "Wild-type",               'T','G','G','T','G','A','T','A'),
        new Mc1rHaplotype("E1",   "Extended black",          'C','A','G','T','G','A','T','A'),
        new Mc1rHaplotype("E2",   "Extended black",          'C','A','G','T','G','A','C','A'),
        new Mc1rHaplotype("E/R1", "Extended black/Birchen",  'T','A','A','T','G','A','C','A'),
        new Mc1rHaplotype("E/R2", "Extended black/Birchen",  'T','A','G','T','G','A','C','A'),
        new Mc1rHaplotype("R1",   "Birchen",                 'T','A','G','T','G','A','T','A'),
        new Mc1rHaplotype("R2",   "Birchen",                 'T','G','G','A','G','A','C','A'),
        new Mc1rHaplotype("R3",   "Birchen",                 'T','G','G','A','G','A','T','A'),
        new Mc1rHaplotype("R4",   "Birchen?",                'T','A','G','C','G','A','C','A'),
        new Mc1rHaplotype("B1",   "Brown",                   'T','A','G','T','A','A','T','A'),
        new Mc1rHaplotype("B2",   "Brown",                   'T','G','G','T','G','A','T','C'),
        new Mc1rHaplotype("B3",   "Brown?",                  'T','A','G','T','G','A','T','C'),
        new Mc1rHaplotype("B/BC", "Brown/Buttercup",         'C','A','G','T','G','A','T','C'),
        new Mc1rHaplotype("BC1",  "Buttercup",               'T','G','G','C','G','G','C','A'),
        new Mc1rHaplotype("BC2",  "Buttercup?",              'T','G','G','C','G','A','C','A'),
        new Mc1rHaplotype("Y1",   "Wheaten",                 'T','G','G','T','G','G','C','A'),
        new Mc1rHaplotype("Y2",   "Wheaten",                 'T','G','G','T','G','G','T','A')
    };

    public static string NormalizeSampleId(string fileName)
    {
        string stem = Path.GetFileNameWithoutExtension(fileName).Trim();
        stem = Regex.Replace(stem, @"[-_.]\d+[-_.](?:P[-_.]?[FR]|EXT(?:ERNAL)?[-_.]?F)$", "", RegexOptions.IgnoreCase);
        stem = Regex.Replace(stem, @"[-_.](?:EXT(?:ERNAL)?[-_.]?F|INT(?:ERNAL)?[-_.]?[FR]|INTERNAL[-_.]?[FR]|[FR])$", "", RegexOptions.IgnoreCase);
        return stem.TrimEnd('-', '_', '.');
    }

    public static SampleConsensusResult Interpret(string sampleId, IReadOnlyList<Mc1rRead2026Result> reads)
    {
        if (reads.Count == 0) return EmptyNoCall(sampleId, "No reads supplied.");
        var usable = reads.Where(r => !r.IsDirty).ToList();
        if (usable.Count == 0) return EmptyNoCall(sampleId, "All reads failed global chromatogram QC.", reads.Count);

        var c212 = BuildConsensus(Mc1r2026Caller.SiteC212, usable.Select(r => r.C212));
        var c274 = BuildConsensus(Mc1r2026Caller.SiteC274, usable.Select(r => r.C274));
        var c376 = BuildConsensus(Mc1r2026Caller.SiteC376, usable.Select(r => r.C376));
        var c398 = BuildConsensus(Mc1r2026Caller.SiteC398, usable.Select(r => r.C398));
        var c409 = BuildConsensus(Mc1r2026Caller.SiteC409, usable.Select(r => r.C409));
        var c427 = BuildConsensus(Mc1r2026Caller.SiteC427, usable.Select(r => r.C427));
        var c637 = BuildConsensus(Mc1r2026Caller.SiteC637, usable.Select(r => r.C637));
        var c644 = BuildConsensus(Mc1r2026Caller.SiteC644, usable.Select(r => r.C644));

        var sites = new[] { c212, c274, c376, c398, c409, c427, c637, c644 };
        bool anyDiscordance = sites.Any(s => s.IsDiscordant);
        int calledSites = sites.Count(s => s.Genotype != "NoCall");
        bool incomplete = calledSites < sites.Length;

        var observed = sites.Where(s => s.Genotype != "NoCall").ToDictionary(s => s.CdnaPosition, s => s.Genotype);
        var candidates = FindCompatibleDiplotypes(observed);

        string callStatus;
        string breedingCategory;
        string interpretation;

        if (candidates.Count == 0)
        {
            callStatus = "UNRESOLVED";
            breedingCategory = "UNRESOLVED";
            interpretation = "Observed MC1R genotype is not represented by any pair of the 18 common Ma et al. (2026) haplotypes. Review chromatograms and consider an uncommon/other haplotype.";
        }
        else
        {
            var categories = candidates.Select(c => BreedingCategoryFor(c.A, c.B)).Distinct().OrderBy(x => x).ToList();
            breedingCategory = categories.Count == 1 ? categories[0] : "AMBIGUOUS";
            callStatus = incomplete ? "INCOMPLETE" : candidates.Count == 1 ? "RESOLVED" : "AMBIGUOUS";
            interpretation = breedingCategory switch
            {
                "E/E-compatible" => "Both compatible MC1R haplotypes are strict Extended Black E1/E2. This supports a line fixed for the tested Extended Black MC1R haplotypes, subject to QC and assay limitations.",
                "E/ALT" => "One compatible haplotype is strict Extended Black E1/E2 and the other is a different common MC1R haplotype. The bird is not genetically fixed E/E at MC1R.",
                "ALT/ALT" => "No compatible diplotype contains a strict Extended Black E1/E2 haplotype. Do not report this bird as E/E-fixed by this assay.",
                _ => "More than one breeder-facing category remains compatible with the observed unphased Sanger genotypes. Do not force an E-status call."
            };
        }

        var qcNotes = new List<string>();
        if (reads.Count != usable.Count) qcNotes.Add($"{reads.Count - usable.Count} read(s) excluded by global QC");
        if (anyDiscordance) qcNotes.Add("one or more sites were discordant between reads");
        foreach (var s in sites.Where(s => !string.IsNullOrWhiteSpace(s.Note))) qcNotes.Add($"c.{s.CdnaPosition}: {s.Note}");

        string qcStatus = qcNotes.Count == 0 ? "OK" : "WARN";
        string diplotypes = candidates.Count == 0 ? "None among 18 common haplotypes" : string.Join("; ", candidates.Select(c => $"{c.A.Id}+{c.B.Id}"));
        string phenotypes = candidates.Count == 0 ? "Unresolved" : string.Join("; ", candidates.Select(c => $"{c.A.Phenotype} + {c.B.Phenotype}").Distinct());

        return new SampleConsensusResult(sampleId, reads.Count, usable.Count, qcStatus, string.Join("; ", qcNotes),
            c212, c274, c376, c398, c409, c427, c637, c644,
            callStatus, breedingCategory, diplotypes, phenotypes, interpretation);
    }

    private static SampleConsensusResult EmptyNoCall(string sampleId, string reason, int totalReads = 0)
    {
        ConsensusSiteCall N(int p) => new(p, "NoCall", 0, false, reason);
        return new SampleConsensusResult(sampleId, totalReads, 0, "FAIL", reason,
            N(212), N(274), N(376), N(398), N(409), N(427), N(637), N(644),
            "NoCall", "NoCall", "NoCall", "NoCall", "No breeding decision should be made from this sample.");
    }

    private static ConsensusSiteCall BuildConsensus(int position, IEnumerable<SiteCall> calls)
    {
        var usable = calls.Where(c => !string.Equals(c.Genotype, "NoCall", StringComparison.OrdinalIgnoreCase)).ToList();
        if (usable.Count == 0) return new ConsensusSiteCall(position, "NoCall", 0, false, "No usable read covers this site.");

        var groups = usable.GroupBy(c => NormalizeGenotype(c.Genotype))
            .Select(g => new { Genotype = g.Key, Count = g.Count(), Notes = g.Where(x => !string.IsNullOrWhiteSpace(x.Notes)).Select(x => x.Notes!).Distinct().ToList() })
            .OrderByDescending(g => g.Count).ThenBy(g => g.Genotype).ToList();

        if (groups.Count == 1)
        {
            string note = groups[0].Notes.Count > 0 ? string.Join(" | ", groups[0].Notes) : "";
            if (usable.Count == 1) note = string.IsNullOrWhiteSpace(note) ? "Single-read support." : "Single-read support. " + note;
            return new ConsensusSiteCall(position, groups[0].Genotype, groups[0].Count, false, note);
        }

        if (groups[0].Count >= 2 && groups[0].Count > groups[1].Count)
            return new ConsensusSiteCall(position, groups[0].Genotype, groups[0].Count, true,
                $"Majority consensus accepted ({groups[0].Genotype} in {groups[0].Count}/{usable.Count} usable reads); discordant read(s) present.");

        return new ConsensusSiteCall(position, "NoCall", usable.Count, true, "Discordant genotype calls between usable reads; no majority consensus.");
    }

    private static string NormalizeGenotype(string genotype)
    {
        var parts = genotype.Split('/', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
        if (parts.Length != 2 || parts.Any(p => p.Length != 1)) return genotype;
        var a = parts.Select(p => char.ToUpperInvariant(p[0])).OrderBy(c => c).ToArray();
        return $"{a[0]}/{a[1]}";
    }

    private sealed record Diplotype(Mc1rHaplotype A, Mc1rHaplotype B);

    private static List<Diplotype> FindCompatibleDiplotypes(IReadOnlyDictionary<int, string> observed)
    {
        var result = new List<Diplotype>();
        for (int i = 0; i < CommonHaplotypes.Count; i++)
            for (int j = i; j < CommonHaplotypes.Count; j++)
            {
                var a = CommonHaplotypes[i]; var b = CommonHaplotypes[j];
                if (observed.All(kv => GenotypeMatches(kv.Key, kv.Value, a, b))) result.Add(new Diplotype(a, b));
            }
        return result;
    }

    private static bool GenotypeMatches(int position, string observedGenotype, Mc1rHaplotype a, Mc1rHaplotype b)
    {
        var expected = new[] { BaseAt(a, position), BaseAt(b, position) };
        Array.Sort(expected);
        return string.Equals($"{expected[0]}/{expected[1]}", NormalizeGenotype(observedGenotype), StringComparison.OrdinalIgnoreCase);
    }

    private static char BaseAt(Mc1rHaplotype h, int position) => position switch
    {
        212 => h.C212, 274 => h.C274, 376 => h.C376, 398 => h.C398,
        409 => h.C409, 427 => h.C427, 637 => h.C637, 644 => h.C644,
        _ => throw new ArgumentOutOfRangeException(nameof(position))
    };

    private static string BreedingCategoryFor(Mc1rHaplotype a, Mc1rHaplotype b)
    {
        bool aE = IsStrictExtendedBlack(a); bool bE = IsStrictExtendedBlack(b);
        if (aE && bE) return "E/E-compatible";
        if (aE || bE) return "E/ALT";
        return "ALT/ALT";
    }

    private static bool IsStrictExtendedBlack(Mc1rHaplotype h) => h.Id is "E1" or "E2";
}
