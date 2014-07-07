# sisg-mcmc-opengl-computer-demos

C code for Eric Anderson's MCMC simulations/visualizations that he uses in the SISG course

To compile this up, I have used autotools, but I haven't successfully gotten it to 
`make dist` or `make distcheck`.  So, I just committed the `./configure` and other
stuff and so you ought to be able to compile all this stuff like as shown below.

Note that you must configure with the `--without-x` option on a Mac.

```sh
# simple commands to get and compile on a mac
git clone https://github.com/eriqande/sisg-mcmc-opengl-computer-demos.git
cd sisg-mcmc-opengl-computer-demos
git submodule init
git submodule update
./configure --without-x  # --without-x is for Macs.  Leave that off for *nix
make
```

Further instructions on how to use the compiled programs (as well as supporting files, etc) are forthcoming.

