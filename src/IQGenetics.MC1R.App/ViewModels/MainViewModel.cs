using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Runtime.CompilerServices;
using Microsoft.Win32;
using IQGenetics.MC1R.Core.MC1R;
using IQGenetics.MC1R.Core.Sanger;

namespace IQGenetics.MC1R.App.ViewModels;

public sealed class MainViewModel : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;

    public string VersionString => "MVP v0.1";

    private string _referencePath = "";
    public string ReferencePath
    {
        get => _referencePath;
        set { _referencePath = value; OnPropertyChanged(); UpdateCanRun(); }
    }

    public ObservableCollection<string> InputFiles { get; } = new();
    public ObservableCollection<ResultRow> Results { get; } = new();

    private string _statusText = "";
    public string StatusText { get => _statusText; set { _statusText = value; OnPropertyChanged(); } }

    private double _progressValue = 0;
    public double ProgressValue { get => _progressValue; set { _progressValue = value; OnPropertyChanged(); } }

    private string _progressLabel = "Idle";
    public string ProgressLabel { get => _progressLabel; set { _progressLabel = value; OnPropertyChanged(); } }

    public string FileCountLabel => $"{InputFiles.Count} file(s) selected";

    private bool _canRun = false;
    public bool CanRun { get => _canRun; private set { _canRun = value; OnPropertyChanged(); } }

    private bool _canExport = false;
    public bool CanExport { get => _canExport; private set { _canExport = value; OnPropertyChanged(); } }

    public RelayCommand SelectReferenceCommand { get; }
    public RelayCommand AddFilesCommand { get; }
    public RelayCommand ClearFilesCommand { get; }
    public RelayCommand RunAnalysisCommand { get; }
    public RelayCommand ExportCsvCommand { get; }

    public MainViewModel()
    {
        SelectReferenceCommand = new RelayCommand(SelectReference);
        AddFilesCommand = new RelayCommand(AddFiles);
        ClearFilesCommand = new RelayCommand(ClearFiles);
        RunAnalysisCommand = new RelayCommand(async () => await RunAnalysisAsync(), () => CanRun);
        ExportCsvCommand = new RelayCommand(ExportCsv, () => CanExport);
    }

    private void SelectReference()
    {
        var dlg = new OpenFileDialog
        {
            Title = "Select MC1R reference (FASTA/FNA/FA/TXT)",
            Filter = "FASTA files (*.fna;*.fasta;*.fa;*.txt)|*.fna;*.fasta;*.fa;*.txt|All files (*.*)|*.*",
            Multiselect = false
        };

        if (dlg.ShowDialog() == true)
        {
            ReferencePath = dlg.FileName;
            StatusText = "Reference selected.";
        }
    }

    private void AddFiles()
    {
        var dlg = new OpenFileDialog
        {
            Title = "Select Sanger chromatograms (.ab1)",
            Filter = "AB1 files (*.ab1)|*.ab1|All files (*.*)|*.*",
            Multiselect = true
        };

        if (dlg.ShowDialog() == true)
        {
            foreach (var f in dlg.FileNames)
            {
                if (!InputFiles.Contains(f))
                    InputFiles.Add(f);
            }
            OnPropertyChanged(nameof(FileCountLabel));
            UpdateCanRun();
        }
    }

    private void ClearFiles()
    {
        InputFiles.Clear();
        Results.Clear();
        CanExport = false;
        StatusText = "Cleared.";
        ProgressValue = 0;
        ProgressLabel = "Idle";
        OnPropertyChanged(nameof(FileCountLabel));
        UpdateCanRun();
    }

    private void UpdateCanRun()
    {
        CanRun = !string.IsNullOrWhiteSpace(ReferencePath) && InputFiles.Count > 0;
        RunAnalysisCommand.RaiseCanExecuteChanged();
        ExportCsvCommand.RaiseCanExecuteChanged();
        OnPropertyChanged(nameof(FileCountLabel));
    }

    private async Task RunAnalysisAsync()
    {
        try
        {
            Results.Clear();
            CanExport = false;
            StatusText = "Loading reference...";
            ProgressLabel = "Preparing...";
            ProgressValue = 0;

            var reference = Mc1rReferenceLoader.Load(ReferencePath);

            int n = InputFiles.Count;
            int done = 0;

            foreach (var file in InputFiles)
            {
                done++;
                ProgressValue = 100.0 * (done - 1) / Math.Max(1, n);
                ProgressLabel = $"Analyzing {done}/{n}...";
                StatusText = Path.GetFileName(file);

                await Task.Run(() =>
                {
                    var ab1 = AbifParser.ReadAb1(file);
                    var result = Mc1rCaller.CallAb1(ab1, reference);

                    App.Current.Dispatcher.Invoke(() =>
                    {
                        Results.Add(new ResultRow
                        {
                            SampleName = result.SampleName,
                            FilePath = result.FilePath,
                            Orientation = result.Orientation.ToString(),
                            AlignmentScore = result.AlignmentScore,
                            DirtyFlag = result.IsDirty ? "DIRTY" : "OK",
                            DirtyReason = result.IsDirty ? result.DirtyReason : "",
                            Genotype274 = result.C274.Genotype,
                            EStatus = result.EStatus,
                            Genotype644 = result.C644.Genotype,
                            SuppressionStatus = result.SuppressionStatus
                        });
                    });
                });
            }

            ProgressValue = 100;
            ProgressLabel = "Done";
            StatusText = $"Completed: {Results.Count} result(s).";
            CanExport = Results.Count > 0;
            ExportCsvCommand.RaiseCanExecuteChanged();
        }
        catch (Exception ex)
        {
            StatusText = "Error: " + ex.Message;
            ProgressLabel = "Error";
        }
    }

    private void ExportCsv()
    {
        var dlg = new SaveFileDialog
        {
            Title = "Export CSV",
            Filter = "CSV (*.csv)|*.csv",
            FileName = "mc1r_results.csv"
        };

        if (dlg.ShowDialog() == true)
        {
            CsvWriter.WriteResults(dlg.FileName, Results);
            StatusText = "Exported: " + dlg.FileName;
        }
    }

    private void OnPropertyChanged([CallerMemberName] string? name = null)
    {
        PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        if (name == nameof(InputFiles))
            OnPropertyChanged(nameof(FileCountLabel));
    }
}
