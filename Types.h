//
// Created by Albert Wong on 3/25/18.
//

#ifndef HULLWHITE_TYPES_H
#define HULLWHITE_TYPES_H

enum FdType {
    Explicit,
    Implicit,
    CrankNicolson,
    SOR
};

enum OptionType {
    Call,
    Put
};

#endif //HULLWHITE_TYPES_H
