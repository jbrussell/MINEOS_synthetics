These are instructions for compiling MINEOS and idagrn on systems that use the ARM-based Apple Silicon chip. They were tested on a Mac Mini using the Apple Silicon M2 Pro. The key to successful compilation is using strictly x86_64 (intel) versions via Rosetta. 

In the future, we will try compiling using native arm64, but for now SAC is only readily available for x86_64.

---
### Step 1: Install gfortran for x86_64 using homebrew

a) Ensure that terminal is running Rosetta. Navigate to the Terminal app in the Finder at "Applications/Utilities". Right click on the Terminal and select Get Info. Check the "Open using Rosetta" check-box. Close the Terminal Info. Now the terminal will install tools with Rosetta translation.

b) Install homebrew for x86_64

    arch -x86_64 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Make sure "which brew" returns: /usr/local/bin/brew. This is the x86_64 version of brew.

If instead it returns /opt/homebrew/bin/brew then you have the arm64 version at the top of your path. You will need to move /usr/local/bin higher in the path using this command:

    export PATH=/usr/local/bin:$PATH

Double check that "which brew" now returns /usr/local/bin/brew.

c) Now install gcc, which includes gfortran

    arch -x86_64 brew install gcc

If installed correctly, "which gfortran" should return: /usr/local/bin/gfortran

---
### Step 2: Install sac for Mac (this will be x86_64). This is required for idagrn6 to generate seismograms. MINEOS will still work without this step.

Request the latest version of SAC for macOS: http://ds.iris.edu/ds/nodes/dmc/forms/sac/

Unzip the directory and put it here: /opt/local/sac

---
### Step 3: Compile Fortran libraries

a) Within this git repository, navigate to ./libgfortran and run: 

    rm ./*.a
    sudo ./makelibs.sh

This will remove existing libraries and compile all new libraries compatible with the same architecture as gfortran.

b) Move back one directory and rename ./bin so you can check at the end that all executables were created successfully.

    cd ..
    mv ./bin ./bin_save

c) Next compile MINEOS and idagrn executables. Run the shell script:

    sudo ./makeall_M2chip.sh

If done correctly, all files should appear in ./bin, and you can delete the saved version:

    rm -r ./bin_save
