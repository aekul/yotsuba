# Yotsuba

Yotsuba is a [Mitsuba](https://github.com/mitsuba-renderer/mitsuba) based renderer that uses [Aether](https://github.com/aekul/aether) to generate samples and compute PDFs.

## Learning About Aether/Yotsuba

The directory `src/integrators/myintegrators/` contains implementations of many standard algorithms using Aether.

A more extensive tutorial is available [here](https://people.csail.mit.edu/lukea/aetherlang).

## Compiling

Yotsuba depends on [Aether](https://github.com/aekul/aether).  Aether is a header-only library, so you can either copy the aether folder to include, or clone the repository directly through the following commands:
```
cd include
git clone https://github.com/aekul/aether
```

Yotsuba uses the same build system as [Mitsuba](https://github.com/mitsuba-renderer/mitsuba) (see Mitsuba's [documentation](https://www.mitsuba-renderer.org/releases/current/documentation.pdf) for build instructions) but there are some things to be aware of when editing `config.py`.

Yotsuba has been tested with clang 3.7 (using -std=c++1z and -std=libc++) so make sure that `config.py` is configured to use an appropriate version of clang and contains the flags:
```bash
-std=c++1z
-std=libc++
```

Yotsuba uses `-std=libc++`, but Mitsuba uses `-std=libstdc++`. Any conflicts between `libc++` and `libstdc++` have been resolved in the source, but have not been resolved in Mitsuba's dependencies. When compiling Yotsuba, you may encounter link errors because of conflicts between `libc++` and `libstdc++`. To resolve these errors, recompile the dependencies manually using `libc++`.

Aether makes extensive use of template metaprogramming. To avoid template/constexpr maximum step/depth related compiler errors, add the following flags to `config.py`:
```bash
-ftemplate-depth=512
-fconstexpr-steps=127124200
-fconstexpr-depth=720
```

## Contact

Please contact Tzu-Mao Li (tzumao at mit.edu) or Luke Anderson (lukea at mit.edu) if there are any issues/comments/questions.
