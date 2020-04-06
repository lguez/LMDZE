# LMDZE
Simplified version of the LMDZ atmospheric general circulation model

## What is it?

LMDZE is a fork of [LMDZ](http://lmdz.lmd.jussieu.fr). The fork dates
from [revision 691 of LMDZ](http://lmdz.lmd.jussieu.fr) in
April 2006. Note that that was the time just before LMDZ was
parallelized. Since then, LMDZE has evolved mostly in form and very
little in content. So, essentially, the content of LMDZE is the
content of LMDZ that was used for
[CMIP3](https://www.wcrp-climate.org/wgcm-cmip/wgcm-cmip3) (or,
equivalently, [AR4](https://www.ipcc.ch/report/ar4/syr)). Actually,
the content is even simpler than the CMIP3 version of LMDZ because the
coupling parts of LMDZ have been removed: LMDZE can only be used with
forced ocean and sea ice, without
[Orchidée](https://orchidee.ipsl.fr). In form, LMDZE has evolved
aiming at clarity and robustness of the source code, taking advantage
of modern features of the Fortran language.

## What is it good for?

What makes the code of an atmospheric general circulation model
("GCM") complex and hard to understand? Many things. The physical
parameterizations and the numerical schemes for solving partial
differential equations may be complex by themselves. If your main
interest is in following and understanding physical parameterizations
and numerical schemes, it is a great help to eliminate other sources
of complexity. Experience shows that parallelism and coupling with
other models bring a code to a whole new level of darkness. Another
constraint that real GCMs have and that significantly weighs on their
clarity is that they must be very careful about backward
compatibility. This means the code becomes loaded with very
little-used options, branchings and alternate versions of
subroutines. So LMDZE breaks free of all this: no parallelism, no
coupling, freedom to eliminate stale code. What remains is a working
GCM that is still close in structure to LMDZ, much simpler to read,
and can thus be used to understand what is going on in LMDZ, or, more
generally, what is in a GCM. Also, LMDZE is a test of ideas on the best
architecture for a large code in modern Fortran.
