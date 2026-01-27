using System;
using System.Collections.Generic;

namespace IQGenetics.MC1R.Core.Alignment;

/// <summary>Simple Smith–Waterman local alignment to map reference positions to query positions.</summary>
public static class SmithWaterman
{
    // Traceback directions
    private const byte STOP = 0;
    private const byte DIAG = 1;
    private const byte UP = 2;
    private const byte LEFT = 3;

    public sealed record Result(int Score, int StartRef, int EndRef, int StartQuery, int EndQuery, Dictionary<int,int> RefToQueryMap);

    public static Result Align(string reference, string query, int match = 2, int mismatch = -1, int gap = -2)
    {
        int n = reference.Length;
        int m = query.Length;

        // score matrix (short is enough for these sizes)
        var H = new short[n + 1, m + 1];
        var TB = new byte[n + 1, m + 1];

        short best = 0;
        int bestI = 0, bestJ = 0;

        for (int i = 1; i <= n; i++)
        {
            char r = reference[i - 1];
            for (int j = 1; j <= m; j++)
            {
                char q = query[j - 1];
                int s = (r == q) ? match : mismatch;

                int diag = H[i - 1, j - 1] + s;
                int up = H[i - 1, j] + gap;
                int left = H[i, j - 1] + gap;

                int val = Math.Max(0, Math.Max(diag, Math.Max(up, left)));
                H[i, j] = (short)val;

                if (val == 0) TB[i, j] = STOP;
                else if (val == diag) TB[i, j] = DIAG;
                else if (val == up) TB[i, j] = UP;
                else TB[i, j] = LEFT;

                if (val > best)
                {
                    best = (short)val;
                    bestI = i;
                    bestJ = j;
                }
            }
        }

        // Traceback to build mapping (only DIAG steps map ref->query)
        var map = new Dictionary<int, int>(capacity: Math.Min(n, m));
        int ti = bestI, tj = bestJ;

        while (TB[ti, tj] != STOP)
        {
            var dir = TB[ti, tj];
            if (dir == DIAG)
            {
                map[ti - 1] = tj - 1;
                ti--; tj--;
            }
            else if (dir == UP)
            {
                ti--;
            }
            else // LEFT
            {
                tj--;
            }
        }

        return new Result(best, ti, bestI, tj, bestJ, map);
    }
}
