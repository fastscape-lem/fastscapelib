/**
 * @file
 * @brief Defines a number of constants used in this library.
 */

#pragma once

namespace fastscapelib
{

    namespace consts
    {

        /*
         * @brief Row index offsets (from a central grid point) of D8
         *  neighbour directions.
         *
         *  The first element correspond to the central grid point itself.
         */
        const short d8_row_offsets[9] = { 0, -1, -1, 0, 1, 1, 1, 0, -1 };


        /**
         * @brief Column index offsets (from a central grid point) of D8
         * neighbour directions.
         *
         * The first element correspond to the central grid point itself.
         */
        const short d8_col_offsets[9] = { 0, 0, -1, -1, -1, 0, 1, 1, 1 };


        /**
         * @brief Row index offsets (from a central grid point) of D8
         *  neighbour directions.
         *
         *  Discard the first element, and only take four neighbors:
         *
         *
         *  \ | /
         *  . X -
         *  . . .
         *
         */
        const short d8_row_offsets_half[4] = { -1, -1, -1, 0 };


        /**
         * @brief Column index offsets (from a central grid point) of D8
         * neighbour directions.
         *
         *  Discard the first element, and only take the first four neighbors
         *
         *
         *  \ | /
         *  . X -
         *  . . .
         *
         */
        const short d8_col_offsets_half[4] = { -1, 0, 1, 1 };

        /**
         * @brief Row index offsets (from a central grid point) of D4
         *  neighbour directions.
         *
         *  The first element correspond to the central grid point itself.
         */
        const short d4_row_offsets[5] = { 0, -1, 0, 0, 1 };


        /**
         * @brief Column index offsets (from a central grid point) of D4
         * neighbour directions.
         *
         * The first element correspond to the central grid point itself.
         */
        const short d4_col_offsets[5] = { 0, 0, -1, 1, 0 };

        /**
         * @brief Row index offsets (from a central grid point) of D4
         *  neighbour directions in cross pattern.
         *
         *  The first element correspond to the central grid point itself.
         */
        const short d4_row_offsets_cross[5] = { 0, -1, -1, 1, 1 };


        /**
         * @brief Column index offsets (from a central grid point) of D4
         * neighbour directions in cross pattern.
         *
         * The first element correspond to the central grid point itself.
         */
        const short d4_col_offsets_cross[5] = { 0, -1, 1, -1, 1 };
    }  // namespace consts

}  // namespace fastscapelib
