# Grid Input Guidelines for Cloudy

When configuring logarithmic density/temperature grids, the parser must see the `hden`/`constant temperature` commands *immediately* before the matching `grid` lines. Cloudy flags “Hydrogen density MUST be specified” if a comment, title, or unrelated directive appears between them.

## Example (works)
```plain
hden -6 vary
grid range from -6 to 6 step 0.5 ncpus 8

constant temperature 4 vary
grid range from 1 to 9 step 0.05
```

## Example (fails)
```plain
hden -6 vary
# log10 n_H from -6 to 6
grid range from -6 to 6 step 0.5 ncpus 8
```

The comment prevents `ParseGrid` from seeing the preceding `vary` command, so hydrogen density stays unset and the run aborts.

### Recommendations
- Keep `grid` lines directly adjacent to their `… vary` commands.
- Move explanatory comments elsewhere (e.g., at the top or end of the file).
- If you need a `title`, place it after the grid section.
- After editing, scan the echoed input in the `.out` file to confirm the parser saw the commands in the intended order.
