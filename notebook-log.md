Recreating the analysis of Glemin et al. (2019) Pervasive hybridizations in the history of wheat relatives.

### Getting the data ###
I navigated to https://www.agap-ge2pop.org/wheat-relative-history/ in the acknowledgements of the paper, downloaded the full archive, and moved it to and unzipped everything in the data/ folder.

### Setting up required software ###
Setting up a software folder
- This will contain all my executables so I will only need to add that folder to my PATH
- Opened up my \~/.bashrc file (\~ is shorthand for /home/ntlin/ folder) -- this configures shell setup things on startup
- Added this line to .bashrc: export PATH="$PATH:~/phylo-practicum/software"
- Now all executables in ~/phylo-practicum/software should be accessible

RAxML
- Downloaded the Linux binary at https://github.com/amkozlov/raxml-ng/releases/download/1.2.2/raxml-ng_v1.2.2_linux_x86_64.zip
- Moved it to software/, unzipped it, and moved raxml-ng from the zip folder to software/
- Did `sudo chmod +x raxml-ng` to make it an executable
- Now I can use the command `raxml-ng`
