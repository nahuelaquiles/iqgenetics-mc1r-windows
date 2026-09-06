using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Runtime.CompilerServices;
using Microsoft.Win32;
using IQGenetics.MC1R.Core.MC1R;
using IQGenetics.MC1R.Core.Sanger;
using System.Threading.Tasks;

namespace IQGenetics.MC1R.App.ViewModels
{
    public sealed class MainViewModel : INotifyPropertyChanged
    {
        public event PropertyChangedEventHandler? PropertyChanged;
        public string VersionString => "v1.0.1 — 2026 haplotype model";

        private string _referencePath = "";
        public string ReferencePath { get => _referencePath; set { _referencePath = value; OnPropertyChanged(); UpdateCanRun(); } }
        public ObservableCollection<string> InputFiles { get; } = new();
        public ObservableCollection<ResultRow> Results { get; } = new();

        private string _statusText = "";
        public string StatusText { get => _statusText; set { _statusText = value; OnPropertyChanged(); } }
        private double _progressValue;
        public double ProgressValue { get => _progressValue; set { _progressValue = value; OnPropertyChanged(); } }
        private string _progressLabel = "Idle";
        public string ProgressLabel { get => _progressLabel; set { _progressLabel = value; OnPropertyChanged(); } }
        public string FileCountLabel => $"{InputFiles.Count} AB1 file(s) selected";
        private bool _canRun;
        public bool CanRun { get => _canRun; private set { _canRun = value; OnPropertyChanged(); } }
        private bool _canExport;
        public bool CanExport { get => _canExport; private set { _canExport = value; OnPropertyChanged(); } }

        public RelayCommand SelectReferenceCommand { get; }
        public RelayCommand AddFilesCommand { get; }
        public RelayCommand ClearFilesCommand { get; }
        public RelayCommand RunAnalysisCommand { get; }
        public RelayCommand ExportCsvCommand { get; }

        public MainViewModel()
        {
            string bundledReference = Path.Combine(AppContext.BaseDirectory, "MC1R_reference.fna");
            if (File.Exists(bundledReference)) _referencePath = bundledReference;
            SelectReferenceCommand = new RelayCommand(SelectReference);
            AddFilesCommand = new RelayCommand(AddFiles);
            ClearFilesCommand = new RelayCommand(ClearFiles);
            RunAnalysisCommand = new RelayCommand(async () => await RunAnalysisAsync(), () => CanRun);
            ExportCsvCommand = new RelayCommand(ExportCsv, () => CanExport);
            UpdateCanRun();
            if (!string.IsNullOrWhiteSpace(_referencePath)) StatusText = "Bundled MC1R reference loaded.";
        }

        private void SelectReference()
        {
            var dlg = new OpenFileDialog { Title = "Select MC1R reference (FASTA/FNA/FA/TXT)", Filter = "FASTA files (*.fna;*.fasta;*.fa;*.txt)|*.fna;*.fasta;*.fa;*.txt|All files (*.*)|*.*", Multiselect = false };
            if (dlg.ShowDialog() == true) { ReferencePath = dlg.FileName; StatusText = "Reference selected."; }
        }

        private void AddFiles()
        {
            var dlg = new OpenFileDialog { Title = "Select all Sanger chromatograms for the batch (.ab1)", Filter = "AB1 files (*.ab1)|*.ab1|All files (*.*)|*.*", Multiselect = true };
            if (dlg.ShowDialog() == true)
            {
                foreach (var f in dlg.FileNames) if (!InputFiles.Contains(f)) InputFiles.Add(f);
                OnPropertyChanged(nameof(FileCountLabel)); UpdateCanRun();
            }
        }

        private void ClearFiles()
        {
            InputFiles.Clear(); Results.Clear(); CanExport = false; StatusText = "Cleared."; ProgressValue = 0; ProgressLabel = "Idle";
            OnPropertyChanged(nameof(FileCountLabel)); UpdateCanRun();
        }

        private void UpdateCanRun()
        {
            CanRun = !string.IsNullOrWhiteSpace(ReferencePath) && InputFiles.Count > 0;
            RunAnalysisCommand?.RaiseCanExecuteChanged(); ExportCsvCommand?.RaiseCanExecuteChanged(); OnPropertyChanged(nameof(FileCountLabel));
        }

        private async Task RunAnalysisAsync()
        {
            try
            {
                Results.Clear(); CanExport = false; StatusText = "Loading reference..."; ProgressLabel = "Preparing..."; ProgressValue = 0;
                var reference = Mc1rReferenceLoader.Load(ReferencePath);
                var readResults = new List<Mc1rRead2026Result>(); int n = InputFiles.Count;
                for (int i = 0; i < n; i++)
                {
                    string file = InputFiles[i]; ProgressValue = 100.0 * i / Math.Max(1, n); ProgressLabel = $"Reading {i + 1}/{n}..."; StatusText = Path.GetFileName(file);
                    var readResult = await Task.Run(() => { var ab1 = AbifParser.ReadAb1(file); return Mc1r2026Caller.CallAb1(ab1, reference); });
                    readResults.Add(readResult);
                }

                ProgressLabel = "Building sample consensus...";
                var groups = readResults.GroupBy(r => Mc1rHaplotypeInterpreter.NormalizeSampleId(r.SampleName), StringComparer.OrdinalIgnoreCase).OrderBy(g => g.Key, StringComparer.OrdinalIgnoreCase).ToList();
                foreach (var group in groups)
                {
                    var result = Mc1rHaplotypeInterpreter.Interpret(group.Key, group.ToList());
                    Results.Add(new ResultRow
                    {
                        SampleName = result.SampleId, TotalReads = result.TotalReads, UsableReads = result.UsableReads, QcStatus = result.QcStatus, QcNotes = result.QcNotes,
                        Genotype212 = result.C212.Genotype, Genotype274 = result.C274.Genotype, Genotype376 = result.C376.Genotype, Genotype398 = result.C398.Genotype,
                        Genotype409 = result.C409.Genotype, Genotype427 = result.C427.Genotype, Genotype637 = result.C637.Genotype, Genotype644 = result.C644.Genotype,
                        CallStatus = result.CallStatus, BreedingCategory = result.BreedingCategory, CompatibleDiplotypes = result.CompatibleDiplotypes,
                        CompatiblePhenotypes = result.CompatiblePhenotypes, Interpretation = result.Interpretation
                    });
                }
                ProgressValue = 100; ProgressLabel = "Done"; StatusText = $"Completed: {Results.Count} bird/sample result(s) from {readResults.Count} chromatogram(s).";
                CanExport = Results.Count > 0; ExportCsvCommand.RaiseCanExecuteChanged();
            }
            catch (Exception ex) { StatusText = "Error: " + ex.Message; ProgressLabel = "Error"; }
        }

        private void ExportCsv()
        {
            var dlg = new SaveFileDialog { Title = "Export sample-level MC1R haplotype CSV", Filter = "CSV (*.csv)|*.csv", FileName = "mc1r_haplotype_results.csv" };
            if (dlg.ShowDialog() == true) { CsvWriter.WriteResults(dlg.FileName, Results); StatusText = "Exported: " + dlg.FileName; }
        }

        private void OnPropertyChanged([CallerMemberName] string? name = null) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
    }
}
