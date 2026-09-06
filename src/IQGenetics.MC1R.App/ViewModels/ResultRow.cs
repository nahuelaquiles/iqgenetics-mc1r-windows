namespace IQGenetics.MC1R.App.ViewModels
{
    public sealed class ResultRow
    {
        public string SampleName { get; set; } = "";
        public int TotalReads { get; set; }
        public int UsableReads { get; set; }
        public string QcStatus { get; set; } = "";
        public string QcNotes { get; set; } = "";
        public string Genotype212 { get; set; } = "";
        public string Genotype274 { get; set; } = "";
        public string Genotype376 { get; set; } = "";
        public string Genotype398 { get; set; } = "";
        public string Genotype409 { get; set; } = "";
        public string Genotype427 { get; set; } = "";
        public string Genotype637 { get; set; } = "";
        public string Genotype644 { get; set; } = "";
        public string CallStatus { get; set; } = "";
        public string BreedingCategory { get; set; } = "";
        public string CompatibleDiplotypes { get; set; } = "";
        public string CompatiblePhenotypes { get; set; } = "";
        public string Interpretation { get; set; } = "";
    }
}
