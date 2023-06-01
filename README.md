# NGS Primer Design

### Sherwood Lab

**06/01/2023**

## Goal

Design primers for use in multiplex sequencing. Ideally, the primer set will not contain primers that easily cross-dimerize with any other primer in the pool.

## Overview
1. **Web:** Generate initial primer set on primer3 website (Rich did this step — his standard settings)
2. **R:** Reformat primer3 into fasta format, with the header target|direction|rank, where target and rank are both integer indices
3. **Web/excel:** Utilize PrimerDimer (at PrimerDimer.com) to identify the free energies of dimer formation for all possible primer pairs
    1. Copy and paste fasta file as input 
    2. Select "multiplex analysis" and "excel summary report"
    3. Download report (and be sure to save as csv)
4. **R:** Create symmetric interaction matrix for all primers, with entries as PrimerDimer delta G calculations for all pairs of primers
    1. The data should represent what would be the contents of a symmetric matrix, but should actually be in a long dataframe format as follows:
    2. For each primer pair, there should be 3 columns corresponding to Primer A (target_id index, direction, and rank — rank is referred to as “pair_id) and 3 for Primer B, plus a dG column for the interaction between the two
        1. Each pair should be represented twice in the dataframe, with each primer in the pair appearing once as Primer A and once as Primer B (hence “symmetric”)
5. **R:** Hard filter to remove overall “bad” primers
    1. Identify all primers with an average interaction dG <= -7
    2. Remove rows containing those primers in either the A or B positions (since the symbolic matrix is symmetric, as described above, there is redundancy in the interactions and it is important to remove both occurrences of each bad interaction/primer)
    3. *Generate list of all valid primers (as opposed to valid pairs) for future reference*
6. **R:** Recursive pairwise filtering
    1. Identify the worst (lowest dG) interaction for each primer (target+direction+rank)
    2. From that filtered worst-by-primer set, identify the worst (lowest dG) interaction for each target+direction pair
    3. If still left with multiple pairs for a given target+direction (i.e. a dG tie), choose one randomly
    4. From this list of only worst interactions, select the primer in each pair with the higher (worse) rank/pair_id
    5. This is the list of new “worst” primers to remove. Drop them entirely from the initial matrix
    6. Repeat this process until the worst interaction left in the dataframe is at or exceeds the threshold value (-7)
7. **R:** Run the above process 300 times
    1. For each run, count the number of targets remaining after the filtering process
    2. Return the run/primer set with the most targets remaining
8. *Now it gets tricky — have to fill in the gaps, since some target+direction pairs are removed entirely by the above process…*
9. **R:** Identify targets missing only one direction of their primer pair, and then all the targets missing primers in both directions
    1. Manually go through the original set of valid primers that remained after the hard filtering step (Step 5) to figure out potential primers that could fill in for these missing targets
    2. Add all of the potential options to the dataset of selected primers
    3. Choose the best ranked primer for each target-direction pair
        1. At this point, have only one primer selected from each target-direction pair — 40 rows, 20 targets in two directions each
10. **R:** Identify all interactions (all dGs) across this new complete but reduced set
    1. Inspect this dataframe manually — identify all pairs that have dG below -8 or so
    2. Inspect the complete selected set, and identify all alternative choices for the primers that were in these bad pairs
        1. Most targets will have multiple options in the “selected_primer_set” dataframe generated after the 300 run-and-pick-best step (Step 7) — which other primers in the same direction are in this set but not in the bad interaction category above?
        2. Try replacing a primer from each bad set with another valid option.
    3. Recurse — identify all interactions across new set, identify pairs with bad interactions, and keep trying to swap out primers one at a time until the set qualitatively seems “good enough”
        1. I recommend saving a few different sets with some different choices and validating/evaluating in OpenPrimeR…
11. **R — OpenPrimeR:** use the R library OpenPrimeR to visualize interaction dG values for each set of primers generated above
    1. Use the plot_constraint method to generate boxplots visualizing all interaction dG terms found within each primer set. 
    2. OpenPrimeR has a default minimum dG of (-6/-7) — the plot_constraint method also shows a horizontal line across this minimum value. To easily check your set, just see if its boxplot is fully above this line.

## Results
* I generated 4 different options, and all passed the OpenPrimeR threshold in the plot_constraints visualization — picked one set at random.
* Replicating: current method that wraps plot_constraints  won’t run…
    * I made the mistake of reformatting Set 1 so the column names don’t align — but it should run if you remove the Set 1 csv from the results directory or reformat it back to original style
    * I also haven’t yet successfully installed OpenPrimeR and all its dependencies on the lab laptop.
        * openPrimeR is available through bioconductor: https://www.bioconductor.org/packages/release/bioc/html/openPrimeR.html 
        * Dependency needed for dimerization energies (OligoArrayAux) is no longer maintained or available for download

## Contributors

<!-- ALL-CONTRIBUTORS-LIST:START -->
| Contributions | Name |
| ----: | :---- |
| [💻](# "Code") [🤔](# "Ideas and Planning") [📖](# "Documentation")| [Lara Brown](https://github.com/larabbrown) |
| [🔣](# "Data") | Jake Francoeur | 
| [📆](# "Project Management") [🤔](# "Ideas and Planning") | Richard Sherwood |
| [🤔](# "Ideas and Planning") | [Chris Cassa](https://github.com/cassalab) |

<!-- ALL-CONTRIBUTORS-LIST:END -->

(For a key to the contribution emoji or more info on this format, check out [“All Contributors.”](https://allcontributors.org/docs/en/emoji-key))

## How to Provide Feedback

Questions, bug reports, and feature requests can be submitted to this repo's [issue queue](https://github.com/TYTYBU/vcfByGene/issues).

## Questions

Please reach out to any of the contributors above with any questions you may have.