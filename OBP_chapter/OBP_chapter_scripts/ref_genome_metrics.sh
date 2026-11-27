# Get summary of ref genome metrics
# Get scaffold lengths and cumulative positions from the reference genome index
awk '
BEGIN {
    OFS="\t"
    print "Scaffold", "Length", "Cumulative_start", "Cumulative_end", \
          "Start_percent", "End_percent"
}

{
    scaffold_length = $2

    cumulative_start = cumulative
    cumulative_end   = cumulative + scaffold_length
    cumulative       = cumulative_end

    names[NR]   = $1
    lengths[NR] = scaffold_length
    starts[NR]  = cumulative_start
    ends[NR]    = cumulative_end
}

END {
    total = cumulative

    for (i = 1; i <= NR; i++) {
        print names[i], \
              lengths[i], \
              starts[i], \
              ends[i], \
              100 * starts[i] / total, \
              100 * ends[i] / total
    }
}
' bNeochry1.pri.20230906.fa.fai > scaffold_reference_positions.tsv

# Later in R
# Install packages if required
# install.packages(c("tidyverse", "scales"))

library(tidyverse)
library(scales)

# Read scaffold positions
scaffolds <- read_tsv(
  "scaffold_reference_positions.tsv",
  show_col_types = FALSE
)

# Preserve the scaffold order from the FASTA index
scaffolds <- scaffolds %>%
  mutate(
    Scaffold = factor(Scaffold, levels = Scaffold),
    group = row_number() %% 2
  )

# Plot the reference genome as one continuous bar
ggplot(scaffolds) +
  geom_rect(
    aes(
      xmin = Cumulative_start,
      xmax = Cumulative_end,
      ymin = 0,
      ymax = 1,
      fill = factor(group)
    ),
    colour = NA
  ) +
  scale_x_continuous(
    labels = label_number(
      scale = 1e-6,
      suffix = " Mb"
    ),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c("grey35", "grey70"),
    guide = "none"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Orange-bellied parrot reference genome",
    subtitle = "Scaffolds arranged in FASTA index order",
    x = "Cumulative genome position",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.title = element_text(face = "bold")
  )



# Calculate length of callable genome
REF="/projects/erode/people/qhk515/reference_genomes/OBP/bNeochry1.pri.20230906.fa.fai"
LIST="OBP_scaffolds_ge100kb.list"

awk '
NR == FNR {
    keep[$1] = 1
    next
}

{
    total += $2

    if ($1 in keep) {
        kept += $2
        n++
    }
}

END {
    printf "Scaffolds in list: %d\n", n
    printf "Length retained: %d bp\n", kept
    printf "Reference length: %d bp\n", total

    if (total > 0)
        printf "Percentage retained: %.2f%%\n", 100 * kept / total
    else
        print "ERROR: reference length is zero; check that REF is a .fai file" > "/dev/stderr"
}
' "$LIST" "$REF"
