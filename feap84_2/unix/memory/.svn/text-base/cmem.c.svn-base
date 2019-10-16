/*
 *-----[--.----+----.----+----.----------------------------------------]
 *      * * F E A P * * A Finite Element Analysis Program
 *
 *      Copyright (c) 1984-2014: Regents of the University of California
 *                               All rights reserved
 *
 *      Purpose: C Language dynamic memory allocation scheme routines
 *
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel (27 Jan 03)
 *                               Revised : David Bindel (23 Jan 06)
 *-----[--.----+----.----+----.----------------------------------------]
 * Modification log                                Date (dd-mm-year)
 *   Original version                                    01-11-2006
 *   1. Add options for underscore/no underscore for
 *      Fortran calling interface options.               27-11-2006
 *   2. Correct realloc (p = q at 277), add line 283     15-01-2009
 *   3. Delete double zero set on memory allocation      03-05-2011
 *      Add line feed to header corruption error message   
 *-----[--.----+----.----+----.----------------------------------------]
 */

       /*********************************/
       /* Memory Allocation Routines    */
       /* Change 7-places for EACH      */
       /* Type of workstation Search    */
       /* for string FORTRAN to change  */
       /* Currently set: Underscore     */
       /*   Version for: Underscore     */
       /*       GCC, INTEL, SUN & DEC   */
       /*   Version for: No Underscore  */
       /*       HP and IBM              */
       /*********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Types for a single byte or a pointer offset. */
#define byte      unsigned char
#define offset_t  long

static byte* mem_system_base_r;  /* Base address of hr (double) array  */
static byte* mem_system_base_i;  /* Base address of mr (integer) array */

unsigned long compute_crc32(byte* buf, int len);


/*
 * Notes:
 *
 * A block allocated using this interface will have the following layout:
 *
 *  +---+--------+---- .... ----+
 *  | x | p      | user space   |
 *  +---+--------+---- .... ----+
 *
 * where x is between zero and eight bytes of padding for alignment
 * purposes; p is a structure containing memory manager data; and
 * "user space" is the memory that the user can access.  The number of
 * bytes in the padding field (x) is chosen so that the beginning of
 * the user space is a multiple of eight bytes distant from the base
 * address of a reference array (hr or mr).
 *
 * The blocks can consist of four-byte integers (atomlen = 4), or
 * eight-byte integers or doubles (atomlen = 8).  The beginning of
 * the user data field for a block is at position
 *
 *     base + atomlen * c_user_offset
 * 
 * In the Fortran code, this offset is at
 *
 *     f_user_offset = c_user_offset + 1
 *
 * The C memory manager structure is at
 *
 *     c_user_offset - sizeof(cmem_data_t)
 *
 * To keep from running into confusing alignment issues when accessing
 * the C memory manager structure, we will always work with a local copy of
 * the user data field.
 *
 * When reallocfn is used to reallocate a block, there may be a different
 * amount of alignment padding needed for the old and the new chunks.
 * Consequently, we look at the old and the new sizes of the padding field x,
 * and shift the contents of the realloc'ed array by the difference in the
 * sizes.
 */


typedef struct {
    void* p_save;
    int   len;
    unsigned long header_checksum;
    unsigned long data_checksum;
} cmem_data_t;

#define cmem_bytes  sizeof(cmem_data_t)
#define cmem_chunks ((cmem_bytes+7)/8)



static unsigned long cmem_header_checksum(cmem_data_t* cmem_data)
/*
 *-------------------------------------------------------
 *      Purpose: Compute a checksum quantity over the header
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    cmem_data_t tmp;
    memcpy(&tmp, cmem_data, cmem_bytes);
    tmp.header_checksum = 0;
    return compute_crc32((byte*) &tmp, cmem_bytes);
}


static void set_cmem_data(byte* base, offset_t offset, int atomlen, 
                          cmem_data_t* cmem_data) 
/*
 *-------------------------------------------------------
 *      Purpose: Set the CMEM data for a block
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    if (offset == -1) {
        printf("Cannot access block at zero offset\n");
        exit(1);
    }
    cmem_data->header_checksum = cmem_header_checksum(cmem_data);
    memcpy(base + offset * atomlen - cmem_bytes, cmem_data, cmem_bytes);
}


static void get_cmem_data(byte* base, offset_t offset, int atomlen,
                          cmem_data_t* cmem_data)
/*
 *-------------------------------------------------------
 *      Purpose: Get the CMEM data for a block
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    if (offset == -1) {
        printf("Cannot access block at zero offset\n");
        exit(1);
    }
    memcpy(cmem_data, base + offset * atomlen - cmem_bytes, cmem_bytes);
    if (cmem_data->header_checksum != cmem_header_checksum(cmem_data)) {
        printf("Header field was corrupted\n");
        exit(1);
    }
}


static void cmem_umark(byte* base, offset_t fortran_offset, int atomlen)
/*
 *-------------------------------------------------------
 *      Purpose: C Compute checksum over a FEAP array
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    offset_t offset = fortran_offset - 1;
    byte*    p      = base + atomlen * offset;
    cmem_data_t cmem_data;

    get_cmem_data(base, offset, atomlen, &cmem_data);
    cmem_data.data_checksum = compute_crc32(p, cmem_data.len * 8);
    set_cmem_data(base, offset, atomlen, &cmem_data);
}


static int cmem_ucheck(byte* base, offset_t fortran_offset, int atomlen)
/*
 *-------------------------------------------------------
 *      Purpose: C Return true if no changes since last checksum
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    offset_t offset = fortran_offset - 1;
    byte*    p      = base + atomlen * offset;
    cmem_data_t cmem_data;

    get_cmem_data(base, offset, atomlen, &cmem_data);
    return (cmem_data.data_checksum == compute_crc32(p, 8*cmem_data.len));
}



static offset_t allocfn(byte* base, int len, int atomlen)
/*
 *-------------------------------------------------------
 *      Purpose: C Dynamic allocate memory for FEAP arrays
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    int atoms_per_8 = 8 / atomlen;

    byte*       p;           /* Start of calloc'ed space */
    offset_t    offset;      /* Offset of user data      */
    cmem_data_t cmem_data;   /* CMEM data field          */

    p       = (byte*) calloc(len + 1 + cmem_chunks, 8);
    offset  = atoms_per_8 * ((offset_t) (p-base + cmem_bytes + 7) / 8);

    if (p == NULL) {
        printf("Memory allocation error\n CALLOC() returns NULL pointer\n");
        return 0;
    } 

    memset(&cmem_data, 0, cmem_bytes);
    cmem_data.p_save        = p;
    cmem_data.len           = len;
    set_cmem_data(base, offset, atomlen, &cmem_data);

    return offset + 1;
}


static void freefn(byte* base, offset_t fortran_offset, int atomlen)
/*
 *-------------------------------------------------------
 *      Purpose: C Free allocated memory for FEAP arrays
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    cmem_data_t cmem_data;
    get_cmem_data(base, fortran_offset-1, atomlen, &cmem_data);
    free(cmem_data.p_save);
}


static offset_t reallocfn(byte* base, offset_t fortran_offset, 
                          int len, int atomlen)
/*
 *-------------------------------------------------------
 *      Purpose: C Dynamic re-allocate memory for FEAP arrays
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *-------------------------------------------------------
 */
{
    int atoms_per_8 = 8 / atomlen;

    byte*       p;            /* Start of old space       */
    offset_t    p_offset;     /* Offset of old user data  */
    cmem_data_t p_cmem_data;  /* Old CMEM data field      */
    int         p_shift;      /* Old user start - base    */

    byte*       q;            /* Start of new space       */
    offset_t    q_offset;     /* Offset of new user data  */
    cmem_data_t q_cmem_data;  /* New CMEM data field      */
    int         q_shift;      /* New user start - base    */

    get_cmem_data(base, fortran_offset-1, atomlen, &p_cmem_data);

    p        = p_cmem_data.p_save;
    p_offset = fortran_offset - 1;
    p_shift  = (base + atomlen * p_offset) - p;

    q        = (byte*) realloc(p, (len + 1 + cmem_chunks)*8);
    q_offset = atoms_per_8 * ((offset_t) (q-base + cmem_bytes + 7) / 8);
    q_shift  = (base + atomlen * q_offset) - q;

    if (q == NULL) {
        printf("Memory allocation error\n CALLOC() returns NULL pointer\n");
        return 0;
    }

    if (p_shift != q_shift)
        memmove(q + q_shift, q + p_shift, len * 8);

    memset(&q_cmem_data, 0, cmem_bytes);
    q_cmem_data.p_save        = q;
    q_cmem_data.len           = len;
    set_cmem_data(base, q_offset, atomlen, &q_cmem_data);

    return q_offset + 1;
}

/* FORTRAN Interface */

/* HP  and IBM use :
void initm(void* base_r, void* base_i)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void initm_(void* base_r, void* base_i)
/*
 *-------------------------------------------------------
 *      Purpose: C Language dynamic memory allocation scheme 
 *               initialization function for FEAP
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         base_r              -- pointer for the base of real arrays
 *         base_i              -- pointer for the base of integer arrays
 *
 *      Output: (stored in local static variables)
 *         mem_system_base_r   -- base pointer for reals
 *         mem_system_base_i   -- base pointer for integers
 *-------------------------------------------------------
 */
{
    mem_system_base_r = (byte*) base_r;
    mem_system_base_i = (byte*) base_i;
}

/* FORTRAN Interface */

/* HP  and IBM use :
void fallocfn(offset_t* offset, int* len, int* precis, int* ipr)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void fallocfn_(offset_t* offset, int* len, int* precis, int* ipr)
/*
 *-------------------------------------------------------
 *      Purpose: Dynamically allocate memory for FEAP arrays
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         len    -- length of desired memory segment in 8-byte words
 *         precis -- (=1) for integers, 
 *                   (=2) for reals
 *         ipr    -- (=1) for 8 byte integers
 *                   (=2) for 4 byte integers
 *
 *      Output: 
 *         offset -- offset to the allocated memory segment
 *                   offset is into hr for precis = 2 and
 *                   into mr for precis = 1
 *-------------------------------------------------------
 */
{
    if (*precis == 1 && *ipr == 2)
        *offset = allocfn(mem_system_base_i, *len, 4);
    else if (*precis == 1 && *ipr == 1)
        *offset = allocfn(mem_system_base_i, *len, 8);
    else if (*precis == 2)
        *offset = allocfn(mem_system_base_r, *len, 8);
    else {
        printf("Memory option precis = %d, ipr = %d not implemented\n",
               *precis, *ipr);
        *offset = 0;
    }
}

/* FORTRAN Interface */

/* HP  and IBM use :
void ffreefn(offset_t* offset, int* precis, int* ipr)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void ffreefn_(offset_t* offset, int* precis, int* ipr)
/*
 *-------------------------------------------------------
 *      Purpose: Free allocated memory for FEAP arrays
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         offset -- offset to memory location to be freed
 *         precis -- (=1) for integers, (=2) for reals
 *         ipr    -- (=1) for 8 byte ints,   (=2) for 4 byte ints
 *-------------------------------------------------------
 */
{
    if (*precis == 1 && *ipr == 2)
        freefn(mem_system_base_i, *offset, 4);
    else if (*precis == 1 && *ipr == 1)
        freefn(mem_system_base_i, *offset, 8);
    else if (*precis == 2)
        freefn(mem_system_base_r, *offset, 8);
    else
        printf("Memory option precis = %d, ipr = %d not implemented\n",
               *precis, *ipr);

    *offset = 0;
}

/* FORTRAN Interface */

/* HP  and IBM use :
void freallocfn(offset_t* offset, int* len, int* precis, int* ipr)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void freallocfn_(offset_t* offset, int* len, int* precis, int* ipr)
/*
 *-------------------------------------------------------
 *      Purpose: Dynamically re-allocate memory for FEAP arrays
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         offset -- current offset
 *         len    -- length of desired memory segment in 8-byte words
 *         precis -- (=1) for integers, 
 *                   (=2) for reals
 *         ipr    -- (=1) for 8 byte integers
 *                   (=2) for 4 byte integers
 *
 *      Output: 
 *         offset -- new offset to the allocated memory segment
 *                   offset is into hr for precis = 2 and
 *                   into mr for precis = 1
 *-------------------------------------------------------
 */
{
    if (*precis == 1 && *ipr == 2)
        *offset = reallocfn(mem_system_base_i, *offset, *len, 4);
    else if (*precis == 1 && *ipr == 1)
        *offset = reallocfn(mem_system_base_i, *offset, *len, 8);
    else if (*precis == 2)
        *offset = reallocfn(mem_system_base_r, *offset, *len, 8);
    else {
        printf("Memory option precis = %d, ipr = %d not implemented\n",
               *precis, *ipr);
        *offset = 0;
    }
}

/* FORTRAN Interface */

/* HP  and IBM use :
void cmemcheck(offset_t* offset, int* precis, int* ipr)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void cmemcheck_(offset_t* offset, int* precis, int* ipr)
/*
 *-------------------------------------------------------
 *      Purpose: Check to ensure memory header is still fine
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         offset -- offset to memory location to be checked
 *         precis -- (=1) for integers, (=2) for reals
 *         ipr    -- (=1) for 8 byte ints,   (=2) for 4 byte ints
 *-------------------------------------------------------
 */
{
    cmem_data_t cmem_data;
    if (*precis == 1 && *ipr == 2)
        get_cmem_data(mem_system_base_i, *offset-1, 4, &cmem_data);
    else if (*precis == 1 && *ipr == 1)
        get_cmem_data(mem_system_base_i, *offset-1, 8, &cmem_data);
    else if (*precis == 2)
        get_cmem_data(mem_system_base_r, *offset-1, 8, &cmem_data);
    else
        printf("Memory option precis = %d, ipr = %d not implemented\n",
               *precis, *ipr);
}

/* FORTRAN Interface */

/* HP  and IBM use :
void cmemumark(offset_t* offset, int* precis, int* ipr)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void cmemumark_(offset_t* offset, int* precis, int* ipr)
/*
 *-------------------------------------------------------
 *      Purpose: Mark the data block with a user checksum
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         offset -- offset to memory location to be checked
 *         precis -- (=1) for integers, (=2) for reals
 *         ipr    -- (=1) for 8 byte ints,   (=2) for 4 byte ints
 *-------------------------------------------------------
 */
{
    if (*precis == 1 && *ipr == 2)
        cmem_umark(mem_system_base_i, *offset, 4);
    else if (*precis == 1 && *ipr == 1)
        cmem_umark(mem_system_base_i, *offset, 8);
    else if (*precis == 2)
        cmem_umark(mem_system_base_r, *offset, 8);
    else
        printf("Memory option precis = %d, ipr = %d not implemented\n",
               *precis, *ipr);

}

/* FORTRAN Interface */

/* HP  and IBM use :
void cmemucheck(offset_t* offset, int* precis, int* ipr, int* result)
*/
/* GCC, INTEL, SUN and DEC use :
*/
void cmemucheck_(offset_t* offset, int* precis, int* ipr, int* result)
/*
 *-------------------------------------------------------
 *      Purpose: Check whether user data has changed since marked
 *               Copyright (c) 2003-2014 : Sanjay Govindjee
 *                                       : David Bindel
 *      Input:
 *         offset -- offset to memory location to be checked
 *         precis -- (=1) for integers, (=2) for reals
 *         ipr    -- (=1) for 8 byte ints,   (=2) for 4 byte ints
 *-------------------------------------------------------
 */
{
    if (*precis == 1 && *ipr == 2)
        *result = cmem_ucheck(mem_system_base_i, *offset, 4);
    else if (*precis == 1 && *ipr == 1)
        *result = cmem_ucheck(mem_system_base_i, *offset, 8);
    else if (*precis == 2)
        *result = cmem_ucheck(mem_system_base_r, *offset, 8);
    else
        printf("Memory option precis = %d, ipr = %d not implemented\n",
               *precis, *ipr);

}


