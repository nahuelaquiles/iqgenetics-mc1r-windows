namespace IQGenetics.MC1R.Core.Sanger;

internal sealed record AbifDirectoryEntry(
    string Name,
    uint Number,
    AbifElementType ElementType,
    ushort ElementSize,
    uint NumElements,
    uint DataSize,
    uint DataOffset,
    uint Handle
);
