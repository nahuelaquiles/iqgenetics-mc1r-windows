using IQGenetics.MC1R.Core.MC1R;

static void Assert(bool condition, string message)
{
    if (!condition) throw new Exception("SELF-TEST FAILED: " + message);
}

static string G(char a, char b)
{
    var v = new[] { a, b }; Array.Sort(v); return $"{v[0]}/{v[1]}";
}

static Mc1rRead2026Result MakeRead(string fileName, Mc1rHaplotype a, Mc1rHaplotype b, string? c644Override = null)
{
    SiteCall S(int p, char x, char y) => new(p, G(x, y), null);
    var qc = new ReadQcSummary(false, "", 700, 600, 30, 0.90, 0.01, 0.01, 0.01, 1, "Clean");
    return new Mc1rRead2026Result(
        fileName, fileName, Orientation.Forward, 1400, false, "",
        S(212, a.C212, b.C212), S(274, a.C274, b.C274), S(376, a.C376, b.C376), S(398, a.C398, b.C398),
        S(409, a.C409, b.C409), S(427, a.C427, b.C427), S(637, a.C637, b.C637),
        new SiteCall(644, c644Override ?? G(a.C644, b.C644), null), qc);
}

static Mc1rHaplotype H(string id) => Mc1rHaplotypeInterpreter.CommonHaplotypes.Single(x => x.Id == id);

var e1 = H("E1"); var e2 = H("E2"); var r1 = H("R1"); var n1 = H("N1"); var n2 = H("N2");

var ee = Mc1rHaplotypeInterpreter.Interpret("EE", new[] { MakeRead("EE-External-F.ab1", e1, e2) });
Assert(ee.BreedingCategory == "E/E-compatible", "E1+E2 must be E/E-compatible.");
Assert(ee.CompatibleDiplotypes.Contains("E1+E2"), "E1+E2 diplotype should be reported.");

var ealt = Mc1rHaplotypeInterpreter.Interpret("EALT", new[] { MakeRead("EALT-Ext-F.ab1", e1, r1) });
Assert(ealt.BreedingCategory == "E/ALT", "E1+R1 must be E/ALT.");

var alt = Mc1rHaplotypeInterpreter.Interpret("ALT", new[] { MakeRead("ALT-External_F.ab1", n1, n2) });
Assert(alt.BreedingCategory == "ALT/ALT", "N1+N2 must be ALT/ALT.");

var majority = Mc1rHaplotypeInterpreter.Interpret("MAJ", new[]
{
    MakeRead("MAJ-1-F.ab1", n1, n1),
    MakeRead("MAJ-1-R.ab1", n1, n1),
    MakeRead("MAJ-1-Ext-F.ab1", n1, n1, "A/C")
});
Assert(majority.C644.Genotype == "A/A", "2:1 consensus must override one discordant c.644 call.");
Assert(majority.QcStatus == "WARN", "Majority consensus with discordance must retain WARN QC.");

// c.212 must never be inferred from Internal-F/Internal-R. This prevents edge artifacts
// from creating biologically impossible c.212 genotypes and forcing a wrong diplotype.
var internalOnly = Mc1rHaplotypeInterpreter.Interpret("INT", new[] { MakeRead("INT-Internal-R.ab1", e1, r1) });
Assert(internalOnly.C212.Genotype == "NoCall", "Internal-R c.212 must be ignored.");
Assert(internalOnly.C212.Note.Contains("External-F"), "Missing External-F should be explicit in c.212 QC note.");
Assert(Mc1rHaplotypeInterpreter.IsExternalForwardReadName("24307-1-2-Ext-F.ab1"), "Ext-F role detection failed.");
Assert(Mc1rHaplotypeInterpreter.IsExternalForwardReadName("24307-1-2-External_F.ab1"), "External_F role detection failed.");
Assert(!Mc1rHaplotypeInterpreter.IsExternalForwardReadName("24307-1-2-P-R.ab1"), "Internal/P-R must not be treated as External-F.");

Assert(Mc1rHaplotypeInterpreter.NormalizeSampleId("24307-1-2-P-F.ab1") == "24307-1", "GENEWIZ P-F naming normalization failed.");
Assert(Mc1rHaplotypeInterpreter.NormalizeSampleId("24307-1-2-Ext-F.ab1") == "24307-1", "External-F naming normalization failed.");
Assert(Mc1rHaplotypeInterpreter.NormalizeSampleId("17663-143JB-B-R.ab1") == "17663-143JB-B", "Simple reverse naming normalization failed.");

Console.WriteLine("MC1R 2026 haplotype self-test: PASS");
