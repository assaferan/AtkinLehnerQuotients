// scripts/automorphisms_645.m
//
// Run directly (from the repo root):  magma -b scripts/automorphisms_645.m
//
// X_0(645)^* is the one genus 5 level with a nontrivial Aut_Q, namely Z/2.
// Checks that its rational points are exactly the special (CM + cusp) points,
// and that the nontrivial automorphism permutes them without producing any new
// point.
//
// phi is built as an iso with its inverse supplied (itself, being an
// involution).  Give Magma the inverse and IsAutomorphism is instant; on a bare
// hom, where it has to find the inverse itself, it does not return in 20
// minutes on this model.

load "src/AtkinLehner.m";           // RationalCMDiscs
load "data/genus5_models.m";

N := 645;
B := 1000;
X := models[N]`curve;
x := X.1; y := X.2; z := X.3; w := X.4; t := X.5;

phi := iso< X -> X | [ z-y, -y, x-y, w-y, x-y-z+t ],
                     [ z-y, -y, x-y, w-y, x-y-z+t ] : Check := true >;

pts := PointSearch(X, B : Nonsingular := true);
cm  := RationalCMDiscs(N);
special := 1 + &+[ cm[d] : d in Keys(cm) ];      // CM points, plus the cusp

printf "X_0(%o)^*: genus %o, degree %o in P^4\n", N, Genus(X), Degree(X);
printf "rational points up to height %o : %o\n", B, #pts;
printf "special points (cusp + CM %o)   : %o\n", Sort([d : d in Keys(cm)]), special;
printf "exceptional points              : %o\n\n", #pts - special;

printf "phi                 : %o\n", DefiningPolynomials(phi);
printf "phi^-1              : %o\n", DefiningPolynomials(Inverse(phi));
printf "IsAutomorphism(phi) : %o\n", IsAutomorphism(phi);
printf "phi ne 1            : %o\n\n", phi ne IdentityMap(X);

for P in pts do
    im := phi(P);
    printf "  %-22o -> %-22o %o\n", Sprint(P), Sprint(im),
        im in pts select "known" else "*** NEW POINT ***";
end for;
printf "\nphi permutes the rational points: %o\n", { phi(P) : P in pts } eq { P : P in pts };
quit;
