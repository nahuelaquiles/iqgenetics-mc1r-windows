namespace IQGenetics.MC1R.App.ViewModels;

public sealed class ResultRow
{
    public string SampleName { get; set; } = "";
    public string FilePath { get; set; } = "";
    public string Orientation { get; set; } = "";
    public int AlignmentScore { get; set; }

    public string DirtyFlag { get; set; } = "";
    public string DirtyReason { get; set; } = "";

    public string Genotype274 { get; set; } = "";
    public string EStatus { get; set; } = "";

    public string Genotype644 { get; set; } = "";
    public string SuppressionStatus { get; set; } = "";
}
