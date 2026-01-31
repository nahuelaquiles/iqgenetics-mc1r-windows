namespace IQGenetics.MC1R.App.ViewModels;

// Representa una fila de resultados en la interfaz
public sealed class ResultRow
{
    public string SampleName { get; set; } = "";
    public string FilePath { get; set; } = "";
    public string Orientation { get; set; } = "";
    public int AlignmentScore { get; set; }
    public string DirtyFlag { get; set; } = "";
    public string DirtyReason { get; set; } = "";
    // Genotipos para cada sitio de interés
    public string Genotype212 { get; set; } = "";
    public string Genotype274 { get; set; } = "";
    public string Genotype355 { get; set; } = "";
    public string Genotype376 { get; set; } = "";
    public string Genotype636 { get; set; } = "";
    public string Genotype637 { get; set; } = "";
    public string Genotype644 { get; set; } = "";
    public string Genotype834 { get; set; } = "";
    public string EStatus { get; set; } = "";
    public string SuppressionStatus { get; set; } = "";
    // Nuevo campo para el patrón de cromatograma segun el QC
    public string ChromatogramPattern { get; set; } = "";
}
