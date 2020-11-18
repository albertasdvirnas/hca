v (future, see CHANGELOG)

v 4.3 
- add bitmask generation for theory
- and it is also used when computing mass pcc
v 4.2 
- mp comparison method is now functional,
- can add consensus

v 4.1 
- add some additional comparison methods (still in development), add posibility to split barcodes into multiple smaller ones
v 4.0.2 (17/03/2020, Current)

v 4.0.1 

v 3.2.8.4
- fix so that barcode generation works when edges only at the borders detected

v 3.2.8.3
- add option to skip edge detection altogether

v 3.2.8.2
-fix bug that wouldn't save large theories

v 3.2.8.1

- fix the bug in the modified otsu code

v 3.2.7

- change the edge detection so that either otsu or tan method would be used, but not both

v 3.2.6

- Made theory generation reproducible by fixing the ring for generating random nucleotides

v 3.2.5

- remove unnecessarily saving of theory files

v 3.2.4 

- fixed up error when plotting barcode vs theory that overlapped to beginning

v 3.2.2

- fix tiny bug with names and when only one timeframe is present
- add two new infoscores

v 3.2.1
v 3.2
- Fix option to add stretching factors
- Output infoscores
- Don't compute comparison for small experiment lengths
- Other minor changes

v 3.1.1.3
- fix the code so it would work with matlab as well

v 3.1.1
- add dialogue to add output folder 
v 3.1.0

- Version with completely reshuffled settings input, and fast generation of results
- In this version it only saves the cc table. 

v 3.0.0
- Added documentation


v 2.1
- reshufled the codes to make the output smaller
- improved the model of stretching to correct nmbp resolution, 
- add function to not remove any timeframes (just select 0 frames)
- add function to output session file after generating barcodes/consensus

v 2.0
- included a faster version for p-values calculation, including possible pre-calculation

v 1.9

v 1.8.2

- save barcode plots

v 1.8.1
- fix bug that prevented output when no filtered barcode selected
- fixed a bug which stretches theory even the nm/bp is the same

v 1.8

- fixed a small indexing typo in the consensus code
- slight update for ssdalign method
- included possibility to change nm/bp ratio after computing theory barcodes
- Some fixes in the visualization of barcode on any chromosome

v 1.7.4

- a number of fixes
- allow inclusion of multiple chromosomes

v 1.6.2

- A number of bug-fixes. No pre-stretching option now works, fixed the bug with wrong bit-masks
- Bitmasks are now always computed before the pre-stretching.


v 1.6

- ssdalign method for alignment now works properly
- fixed small errors in additional results tabs
- changed theory generation into recently updated theory
- Fully implemented two methods for p-value calculation, with first one using extreme value distribution
, second one using the cross correlation coefficient distribution.
The methods work by randomizing the shorter barcode.
- Changed the consensus generation code
- Speed up in comparing theory to experiments as faster implementation of match score comparison is used here.

v 1.5

- added p-value calculation (current version is more or less accurate only for barcodes without stretching, have to update
 the theory for the stretched barcodes.., also the database mean used here is old and should be recalculated, so should add an option to recalculate database average..)
- Added (as an extra results tab) theory barcode vs chromosome comparison, just have the theory barcode generated before.
- Added possibility of using new alignment algorithm ssdalign (parameters for this still need to be fine tuned so it performs worse than 
current method (more alignment errors), but much faster.
- Results tab appears after generating results now
- Added an option to pre-stretch barcodes to the same length, this will make p-value computation faster and increase the number of matching barcodes
- Added an option of filtering before 3% stretching or after for the filtered barcodes. 
- Parallelized kymo alignment and comparison of exp to theory
- rewrote portions of code and added HCA_Code, which can be used as a template for fast comparisons for the same data with small parameters differing
- Moved the consensus barcode of filtered version to the filtered barcode consensus structure
- Added timers to show how much time each step takes

v 1.4.3
- Bugfix

v 1.4.2

- Bugfix

v 1.4.1

- Bugfix

v 1.4

- Filter any number of timeframes, still filtered during the generation
- Fixed consensus result plots, changed markers

v 1.3

- Included different markers '*','x','o' for plots
- Added possibility to change Netropsin/Yoyo concentration when generating theoretical barcode
- hcaSessionStructure is now accessible in the workspace after plotting the results
- Added additional results tab to plot for barcode vs exp plots, p-values, etc.

