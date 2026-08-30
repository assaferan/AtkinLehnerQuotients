// scripts/automorphisms_370.m
//
// Run directly (from the repo root):  magma -b scripts/automorphisms_370.m
//
// X_0(370)^* is the one genus 4 level in this batch, and one of the four
// genus-4 levels that carry an exceptional point.  Its Aut_Q is Z/2.  This
// checks that the nontrivial automorphism permutes the ten rational points
// without producing a new one, and reports where it sends the exceptional
// point in particular.
//
// phi is built as an iso with its inverse supplied (itself, being an
// involution).  Give Magma the inverse and IsAutomorphism is instant; on a bare
// hom, where it has to find the inverse itself, it does not return in 20
// minutes on this model.

load "src/AtkinLehner.m";           // count_special_points_X0Nstar
load "data/genus4_models.m";

N := 370;
B := 1000;
X := models[N]`curve;
x := X.1; y := X.2; z := X.3; w := X.4;

phi := iso< X -> X | [ w-y, -y, -z, x-y ],
                     [ w-y, -y, -z, x-y ] : Check := true >;

// Which rational point is which, from LabelAndAnalyze (src/labeling.m) run on
// the level-370 star forms and carried over to this model by the saturated-HNF
// change of basis.  Keyed by coordinates, so it does not depend on the order
// PointSearch happens to return.
labels := [* <X ! [-1, 1, -1, 1],    "CM -280">,
             <X ! [-1, 2,  0, 1],    "CM -340">,
             <X ! [ 0, 0,  0, 1],    "CM -120">,
             <X ! [ 1, 2,  0, 3],    "CM -36">,
             <X ! [ 0, 0, -1, 1],    "CM -16">,
             <X ! [ 1, 0,  0, 0],    "cusp">,
             <X ! [25,10,  3, 2],    "EXCEPTIONAL">,
             <X ! [ 0, 1, -1, 2],    "CM -40">,
             <X ! [-8,-10,-3,15],    "CM -4">,
             <X ! [ 1, 0,  1, 0],    "CM -100"> *];
exc_pt := X ! [25, 10, 3, 2];

function Label(P)
    for e in labels do
        if e[1] eq P then return e[2]; end if;
    end for;
    return "unlabelled";
end function;

pts := PointSearch(X, B : Nonsingular := true);
special, cm := count_special_points_X0Nstar(N);

printf "X_0(%o)^*: genus %o, degree %o in P^3\n\n", N, Genus(X), Degree(X);
printf "defining equations:\n";
for f in DefiningPolynomials(X) do printf "  %o\n", f; end for;

printf "\nrational points up to height %o : %o\n", B, #pts;
for P in pts do printf "  %-26o %o\n", Sprint(P), Label(P); end for;
printf "special points (cusp + CM %o) : %o\n", Sort([d : d in Keys(cm)]), special;
printf "exceptional points  : %o\n\n", #pts - special;

printf "phi                 : %o\n", DefiningPolynomials(phi);
printf "phi^-1              : %o\n", DefiningPolynomials(Inverse(phi));
printf "IsAutomorphism(phi) : %o\n", IsAutomorphism(phi);
printf "phi ne 1            : %o\n", phi ne IdentityMap(X);
printf "phi^2 = 1           : %o\n\n", forall{ P : P in pts | phi(phi(P)) eq P };

printf "images of the rational points under phi:\n";
for P in pts do
    im := phi(P);
    printf "  %-26o -> %-26o %-14o -> %-14o %o\n", Sprint(P), Sprint(im),
        Label(P), Label(im), im in pts select "known" else "*** NEW POINT ***";
end for;

permutes := { phi(P) : P in pts } eq { P : P in pts };
if permutes then
    G := sub< Sym(#pts) | Sym(#pts) ! [ Index(pts, phi(P)) : P in pts ] >;
    printf "\n|Aut_Q(X_0(%o)^*)| = %o, structure %o\n", N, #G, GroupName(G);
else
    printf "\nphi does not permute the found points: Aut_Q order not read off here\n";
end if;

printf "\nthe exceptional point:\n";
printf "  %o  (%o)\n", exc_pt, Label(exc_pt);
printf "  phi sends it to %o  (%o)\n", phi(exc_pt), Label(phi(exc_pt));
if phi(exc_pt) eq exc_pt then
    printf "  -> fixed by phi\n";
elif phi(exc_pt) in pts then
    printf "  -> swapped with a known rational point, not fixed\n";
else
    printf "  -> *** SENT TO A NEW RATIONAL POINT ***\n";
end if;
printf "  phi does not fix the cusp: cusp -> %o (%o)\n",
    phi(X ! [1,0,0,0]), Label(phi(X ! [1,0,0,0]));

printf "\nphi permutes the rational points: %o\n", permutes;
quit;
