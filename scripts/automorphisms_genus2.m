// scripts/automorphisms_genus2.m
//
// Run directly (from the repo root):  magma -b scripts/automorphisms_genus2.m
//
// Computes Aut_Q(X_0(N)^*) for the seven squarefree genus 2 levels, and checks
// for each whether those automorphisms permute the rational points found up to
// height B or send one of them somewhere the search did not reach.  Prints a
// block per level and a summary table at the end.
//
// Genus 2 levels have no committed model in data/, so the curve comes from
// point_search_X0Nstar, which takes the hyperelliptic branch itself.  Levels
// with no data/starmodels/starforms_<N>.m entry pay for a full Atkin-Lehner
// diagonalization on the first run and get an entry written.

load "src/AtkinLehner.m";

levels := [106, 122, 129, 158, 166, 215, 390];
B := 1000;

rows := [];

for N in levels do
    t0 := Cputime();
    printf "\n========================================================\n";
    printf "N = %o\n", N;

    pts, C := point_search_X0Nstar(N, B);
    Amb<X, Y, Z> := Ambient(C);
    special := count_special_points_X0Nstar(N);

    printf "genus %o, model branch: %o\n", Genus(C),
        IsHyperellipticX0Nstar(N) select "hyperelliptic" else "canonical";
    printf "curve                           : %o\n", C;
    printf "rational points up to height %o : %o\n", B, #pts;
    for P in pts do printf "    %o\n", P; end for;
    printf "special points (cusp + CM)      : %o\n", special;
    printf "exceptional points              : %o\n", #pts - special;

    A, m := AutomorphismGroup(C);
    hyp := HyperellipticInvolution(C);
    printf "\n|Aut_Q| = %o, structure %o (SmallGroup %o)\n", #A, GroupName(A), IdentifyGroup(A);

    newpts := {};
    for a in [g : g in Generators(A) | g ne Id(A)] do
        f := m(a);
        printf "\ngenerator %o : %o%o\n", a, DefiningPolynomials(f),
            f eq hyp select "   [hyperelliptic involution]" else "";
        for P in pts do
            im := f(P);
            printf "  %-22o -> %-22o %o\n", Sprint(P), Sprint(im),
                im in pts select "known" else "*** NEW POINT ***";
            if im notin pts then Include(~newpts, im); end if;
        end for;
    end for;

    // Verdict over the whole group, not just the printed generators.
    permutes := forall{ a : a in A | { m(a)(P) : P in pts } eq { P : P in pts } };
    printf "\nAut_Q permutes the known rational points: %o\n", permutes;
    if #newpts gt 0 then
        printf "*** %o NEW RATIONAL POINT(S) on X_0(%o)^* missed by the height-%o search: %o\n",
            #newpts, N, B, newpts;
    end if;
    printf "[N = %o done in %o s]\n", N, Cputime(t0);

    Append(~rows, <N, #pts, special, #pts - special, #A, GroupName(A), permutes, #newpts>);
end for;

printf "\n========================================================\n";
printf "SUMMARY (genus 2, B = %o)\n\n", B;
printf "%5o %6o %8o %6o %7o %10o %10o %6o\n",
    "N", "#pts", "#special", "#exc", "|Aut|", "Aut", "permutes", "#new";
for r in rows do
    printf "%5o %6o %8o %6o %7o %10o %10o %6o\n",
        r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8];
end for;
printf "\ntotal new rational points found: %o\n", &+[r[8] : r in rows];
quit;
