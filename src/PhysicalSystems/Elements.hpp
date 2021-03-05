#pragma once

enum Element {
    XX,     // 0  = Nothing/vacancy
    H,      // 1  = Hydrogen
    He,     // 2  = Helium
    Li,     // 3  = Lithium
    Be,     // 4  = Beryllium
    B,      // 5  = Boron
    C,      // 6  = Carbon
    N,      // 7  = Nitrogen
    O,      // 8  = Oxygen
    F,      // 9  = Fluorine
    Ne,     // 10 = Neon
    Na,     // 11 = Sodium
    Mg,     // 12 = Magnesium
    Al,     // 13 = Aluminum
    Si,     // 14 = Silicon
    P,      // 15 = Phosphorus
    S,      // 16 = Sulphur
    Cl,     // 17 = Chlorine
    Ar,     // 18 = Argon
    K,      // 19 = Potassium
    Ca,     // 20 = Calcium
    Sc,     // 21 = Scandium
    Ti,     // 22 = Titanium
    V,      // 23 = Vanadium
    Cr,     // 24 = Chromium
    Mn,     // 25 = Manganese
    Fe,     // 26 = Iron
    Co,     // 27 = Cobalt
    Ni,     // 28 = Nickel
    Cu,     // 29 = Copper
    Zn,     // 30 = Zinc
    Ga,     // 31 = Gallium
    Ge,     // 32 = Germanium
    As,     // 33 = Arsenic
    Se,     // 34 = Selenium
    Br,     // 35 = Bromine
    Kr,     // 36 = Krypton
    Rb,     // 37 = Rubidium
    Sr,     // 38 = Strontium
    Y,      // 39 = Yttrium
    Zr,     // 40 = Zirconium
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    I,
    Xe,
    Cs,
    Ba,
    La, 
    Ce, 
    Pr, 
    Nd, 
    Pm, 
    Sm, 
    Eu, 
    Gd, 
    Tb, 
    Dy, 
    Ho, 
    Er, 
    Tm, 
    Yb, 
    Lu,
    Hf, 
    Ta, 
    W, 
    Re, 
    Os, 
    Ir, 
    Pt, 
    Au, 
    Hg, 
    Tl, 
    Pb, 
    Bi, 
    Po, 
    At, 
    Rn,
    Fr, 
    Ra,
    Ac, 
    Th, 
    Pa, 
    U, 
    Np, 
    Pu, 
    Am, 
    Cm, 
    Bk, 
    Cf, 
    Es, 
    Fm, 
    Md, 
    No, 
    Lr,
    Rf, 
    Db, 
    Sg, 
    Bh, 
    Hs, 
    Mt, 
    Ds, 
    Rg, 
    Cn, 
    Nh, 
    Fl, 
    Mc, 
    Lv, 
    Ts, 
    Og
};


Element convertStringToElement(const std::string& str);

std::string convertElementToString(Element& element);
