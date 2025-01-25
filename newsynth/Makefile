SHELL := /bin/bash
GHC_VERSION := 9.4.4
CABAL_VERSION := 3.8.1.0

installghc:
	ghcup install ghc $(GHC_VERSION) CC=g++
	ghcup set ghc $(GHC_VERSION)
	ghcup install cabal $(CABAL_VERSION) CC=g++
	ghcup set cabal $(CABAL_VERSION)
	curl --proto '=https' --tlsv1.2 -sSf https://get-ghcup.haskell.org | sh

buildhs:
	cabal new-build
	cp dist-newstyle/build/aarch64-osx/ghc-$(GHC_VERSION)/newsynth-0.4.0.0/build/gridsynth_ccals/libgridsynth_ccals.dylib .
	cp dist-newstyle/build/aarch64-osx/ghc-$(GHC_VERSION)/newsynth-0.4.0.0/build/libHSnewsynth-0.4.0.0-inplace-ghc$(GHC_VERSION).dylib .

buildtest: buildhs test.c libgridsynth_ccals.dylib libHSnewsynth-0.4.0.0-inplace-ghc$(GHC_VERSION).dylib
	gcc -o test test.c -fPIC -g \
		-I`ghc --print-libdir`/aarch64-osx-ghc-$(GHC_VERSION)/rts-1.0.2/include \
		-L. -lgridsynth_ccals -lHSnewsynth-0.4.0.0-inplace-ghc$(GHC_VERSION) \
		-L`ghc --print-libdir`/aarch64-osx-ghc-$(GHC_VERSION)/ -lHSrts-1.0.2_debug-ghc$(GHC_VERSION)
	
buildtest-with-ghc:
	ghc -o test test.o -no-hs-main -g2 -v\
		-lgridsynth_ccals  -lHSnewsynth-0.4.0.0-inplace-ghc$(GHC_VERSION) -L.

runtest: test
	LD_LIBRARY_PATH=`pwd`:`ghc --print-libdir`/aarch64-osx-ghc-$(GHC_VERSION)/ ./test

