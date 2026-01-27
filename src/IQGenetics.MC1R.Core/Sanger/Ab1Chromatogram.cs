namespace IQGenetics.MC1R.Core.Sanger;

public sealed class Ab1Chromatogram
{
    public string FilePath { get; }
    public string FileName => Path.GetFileName(FilePath);

    public string Bases { get; }
    public int[] Qualities { get; }
    public short[] PeakLocations { get; }
    public string ChannelOrder { get; } // FWO_1, typically "GATC"

    /// <summary>Trace channels keyed by base letter in ChannelOrder. Example: traces['A'] gives the A channel.</summary>
    public IReadOnlyDictionary<char, short[]> Traces { get; }

    public Ab1Chromatogram(string filePath, string bases, int[] qualities, short[] peakLocations, string channelOrder, IReadOnlyDictionary<char, short[]> traces)
    {
        FilePath = filePath;
        Bases = bases;
        Qualities = qualities;
        PeakLocations = peakLocations;
        ChannelOrder = channelOrder;
        Traces = traces;
    }
}
