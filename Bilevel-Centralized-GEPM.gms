$Title Deterministic Single-Node Static Generation Expansion Planning Model (GEPM)
$OnText
                    GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

 Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The GNU General Public License is a free, copyleft license for
software and other kinds of works.

  The licenses for most software and other practical works are designed
to take away your freedom to share and change the works.  By contrast,
the GNU General Public License is intended to guarantee your freedom to
share and change all versions of a program--to make sure it remains free
software for all its users.  We, the Free Software Foundation, use the
GNU General Public License for most of our software; it applies also to
any other work released this way by its authors.  You can apply it to
your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Our General Public Licenses are designed to make sure that you
have the freedom to distribute copies of free software (and charge for
them if you wish), that you receive source code or can get it if you
want it, that you can change the software or use pieces of it in new
free programs, and that you know you can do these things.

  To protect your rights, we need to prevent others from denying you
these rights or asking you to surrender the rights.  Therefore, you have
certain responsibilities if you distribute copies of the software, or if
you modify it: responsibilities to respect the freedom of others.

  For example, if you distribute copies of such a program, whether
gratis or for a fee, you must pass on to the recipients the same
freedoms that you received.  You must make sure that they, too, receive
or can get the source code.  And you must show them these terms so they
know their rights.

  Developers that use the GNU GPL protect your rights with two steps:
(1) assert copyright on the software, and (2) offer you this License
giving you legal permission to copy, distribute and/or modify it.

  For the developers' and authors' protection, the GPL clearly explains
that there is no warranty for this free software.  For both users' and
authors' sake, the GPL requires that modified versions be marked as
changed, so that their problems will not be attributed erroneously to
authors of previous versions.

  Some devices are designed to deny users access to install or run
modified versions of the software inside them, although the manufacturer
can do so.  This is fundamentally incompatible with the aim of
protecting users' freedom to change the software.  The systematic
pattern of such abuse occurs in the area of products for individuals to
use, which is precisely where it is most unacceptable.  Therefore, we
have designed this version of the GPL to prohibit the practice for those
products.  If such problems arise substantially in other domains, we
stand ready to extend this provision to those domains in future versions
of the GPL, as needed to protect the freedom of users.

  Finally, every program is threatened constantly by software patents.
States should not allow patents to restrict development and use of
software on general-purpose computers, but in those that do, we wish to
avoid the special danger that patents applied to a free program could
make it effectively proprietary.  To prevent this, the GPL assures that
patents cannot be used to render the program non-free.

Developed by

   Diego Alejandro Tejada Arango
   dtejada@comillas.edu

   November 23, 2020

$OffText

$OnText
The Bilevel formulation in this code is based on Chapter 3 of the Book:
Investment in Electricity Generation and Transmission
DOI: 10.1007/978-3-319-29501-5
Authors:
 Antonio J. Conejo
 Luis Baringo Morales
 Kazempour Jalal
 Afzal S. Siddiqui
$OffText

sets
 c "candidate generating units" /gc1/
 g "existing  generating units" /ge1/
 d "demands                   " /d1 /
 o "operating conditions      " /o1,o2/
;

parameters
CC   (c) "Production                        cost of candidate generating unit c [$/MWh]" /gc1    25/
CE   (g) "Production                        cost of existing  generating unit g [$/MWh]" /ge1    35/ 
IC   (c) "Annualized investment             cost of candidate generating unit c [$/MW ]" /gc1 70000/
PCmax(c) "Maximum production capacity investment of candidate generating unit c [MW]   " /gc1   500/ 
PEmax(g) "Production capacity                    of existing  generating unit g [MW]   " /ge1   400/

Rho  (o) "Weight of operating condition o [h]"
         /o1 6000
          o2 2760/

price(o) "electricity prices [$/MWh]"
;
table PD (d,o) "load demand d at operating condition o [MW]"
        o1  o2
    d1  290 550
;

variable
 of       "objective function [$]"
 lambda(o) "dual variable of balance equation"

positive variables
 pc  (c,o) "Power produced by candidate generating unit c [MW]"
 pe  (g,o) "Power produced by existing  generating unit g [MW]"
 cmax(c  ) "Capacity       of candidate generating unit c [MW]"
 muEm(g,o) "dual variable of max production of candidate unit c"
 muCm(c,o) "dual variable of max production of existing  unit g"
;

equations
 eTotalCost          "Total cost (operating + investment cost)"
 eBalance      (  o) "balance equation at each operating condition o"
 eCandProdLimit(c,o) "candidate production limit"
 eDualConstE   (g,o) "Dual constraint on existing  generating unit g"
 eDualConstC   (c,o) "Dual constraint on candidate generating unit c"
 eStrongDualEq (  o) "strong duality equalities"
;

eTotalCost.. of =e= sum[o,Rho(o)*[sum[g,CE(g)*pe(g,o)]+sum[c,CC(c)*pc(c,o)]]] + sum[c,IC(c)*cmax(c)] ;

eBalance(o).. sum[g,pe(g,o)]+sum[c,pc(c,o)] =e= sum[d,PD(d,o)] ;

eCandProdLimit(c,o).. pc(c,o) =l= cmax(c) ;

eDualConstE(g,o).. CE(g) - lambda(o) + muEm(g,o) =g= 0 ;
eDualConstC(c,o).. CC(c) - lambda(o) + muCm(c,o) =g= 0 ;

eStrongDualEq(o).. sum[g,CE(g)*pe(g,o)]+sum[c,CC(c)*pc(c,o)] =e= lambda(o)*sum[d,PD(d,o)]-sum[g,muEm(g,o)*PEmax(g)]-sum[c,muCm(c,o)*cmax(c)];

* constraints as bounds
cmax.up(c  ) = PCmax(c);
pe.up  (g,o) = PEmax(g);

Model bilevel_nlp /eTotalCost, eBalance, eCandProdLimit, eDualConstE, eDualConstC, eStrongDualEq /;
Model sglevel_lp  /eTotalCost, eBalance, eCandProdLimit                                          /;

* solve bilevel problem
solve bilevel_nlp using nlp minimizing of;

display pc.l, pe.l, cmax.l, lambda.l;

* solve single level problem
solve sglevel_lp using lp minimizing of;
price(o) = eBalance.m(o) / Rho(o) ;
display pc.l, pe.l, cmax.l, price;

* finding prices without dual information of investment
cmax.fx(c) = cmax.l(c);
solve sglevel_lp using lp minimizing of;
price(o) = eBalance.m(o) / Rho(o) ;
display price;
