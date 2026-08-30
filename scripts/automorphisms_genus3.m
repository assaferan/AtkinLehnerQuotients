// scripts/automorphisms_genus3.m
//
// Run directly (from the repo root):  magma -b scripts/automorphisms_genus3.m
//
// Computes Aut_Q(X_0(N)^*) for the eleven genus 3 levels listed below and
// checks whether those automorphisms permute the rational points found by a
// height-1000 search or send one of them to a point the search missed.
// Every one of these X_0(N)^* is a smooth plane quartic, so
// AutomorphismsOfTernaryQuartic gives the automorphisms as 3x3 matrices over Q;
// with geometric := false (the default) it returns the Q-rational ones only.

load "src/AtkinLehner.m";           // count_special_points_X0Nstar
load "data/genus3_models.m";        // models, and the shared ambient P<x,y,z>

levels := [178, 183, 246, 249, 258, 290, 303, 318, 430, 455, 510];
B := 1000;

// Act by a 3x3 matrix on a point of the plane quartic.
function ApplyMatrix(X, M, Pt)
    v := Eltseq(Pt);
    return X ! [ &+[ M[i][j] * v[j] : j in [1..3] ] : i in [1..3] ];
end function;

summary := [* *];

for N in levels do
    X := models[N]`curve;
    f := DefiningPolynomial(X);
    printf "================ N = %o ================\n", N;
    printf "genus %o plane quartic: %o\n", Genus(X), f;

    pts := PointSearch(X, B : Nonsingular := true);
    special := count_special_points_X0Nstar(N);
    printf "rational points up to height %o : %o\n", B, #pts;
    for Pt in pts do printf "    %o\n", Pt; end for;
    printf "special points (cusp + CM)      : %o", special;
    printf "  [committed value %o]\n", models[N]`special_points;
    printf "exceptional points              : %o\n", #pts - special;

    autos := AutomorphismsOfTernaryQuartic(f);
    G, phi := AutomorphismGroupOfTernaryQuartic(f, autos : explicit := true);
    geom := AutomorphismsOfTernaryQuartic(f : geometric := true);
    Ggeom := AutomorphismGroupOfTernaryQuartic(f, geom);
    printf "|Aut_Q| = %o, structure %o   (geometrically: %o, %o)\n",
        #autos, GroupName(G), #geom, GroupName(Ggeom);

    newpts := {};
    gens := [ phi(G.i) : i in [1..NumberOfGenerators(G)] ];
    for M in gens do
        if M eq Parent(M)!1 then continue; end if;
        printf "\ngenerator\n%o\n", M;
        for Pt in pts do
            im := ApplyMatrix(X, M, Pt);
            known := im in pts;
            if not known then Include(~newpts, im); end if;
            printf "  %-24o -> %-24o %o\n", Sprint(Pt), Sprint(im),
                known select "known" else "*** NEW POINT ***";
        end for;
    end for;

    permutes := &and[ { ApplyMatrix(X, M, Pt) : Pt in pts } eq { Pt : Pt in pts }
                      : M in autos ];
    printf "\nAut_Q permutes the known rational points: %o\n", permutes;
    if #newpts gt 0 then
        printf "*** %o NEW RATIONAL POINT(S) ON X_0(%o)^* ***\n", #newpts, N;
        for Pt in newpts do printf "    %o\n", Pt; end for;
    end if;
    printf "\n";

    Append(~summary, <N, #pts, special, #pts - special, #autos, GroupName(G),
                      GroupName(Ggeom), permutes, #newpts>);
end for;

printf "\n================ summary ================\n";
printf "%-6o %-6o %-8o %-6o %-8o %-8o %-10o %-9o %o\n",
    "N", "#pts", "special", "exc", "|Aut_Q|", "Aut_Q", "Aut_Qbar", "permutes", "new";
for s in summary do
    printf "%-6o %-6o %-8o %-6o %-8o %-8o %-10o %-9o %o\n",
        s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9];
end for;
quit;
