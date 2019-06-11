# sisg-mcmc-opengl-computer-demos

C code for Eric Anderson's MCMC simulations/visualizations that he uses in the SISG course

To compile this up, I have used autotools, but I haven't successfully gotten it to 
`make dist` or `make distcheck`.  So, I just committed the `./configure` and other
stuff and so you ought to be able to compile all this stuff like as shown below.

I'm a total noob with autotools.  If anyone know how to do this correctly, then by all
means please fork it, fix it and send me a pull request.

Note that you must configure with the `--without-x` option on a Mac.

The setup below creates a directory called ~/Library/gfmcmc where the programs look for the
GFMCMC "views" files.

```sh
# simple commands to get and compile on a mac
git clone https://github.com/eriqande/sisg-mcmc-opengl-computer-demos.git
cd sisg-mcmc-opengl-computer-demos
git submodule init
git submodule update
mkdir -p ~/Library/gfmcmc
cp ./MCMC_Demos/* ~/Library/gfmcmc/
./configure --bindir=$(pwd)/MCMC_Demos --without-x  # --without-x is for Macs.  Leave that off for Linux
make
make install
```

That will put all the program binaries into the directory MCMC_Demos.  This then should be a directory with all the
programs and data files that are needed to run the demos.  The Views files necessary will have been put into ~/Library/gfmcmc


Further instructions on how to use the compiled programs (as well as supporting files, etc) are forthcoming.

Note that I have gotten this to compile on macs pretty easily.  But haven't yet successfully compiled it on
Linux.
