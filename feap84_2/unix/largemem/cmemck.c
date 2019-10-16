/*
 * Checksum code from RFC 1952 (GZIP format specification)
 *----[--.----+----.----+----.-----------------------------------------]
 *    Modification log                                Date (dd-mm-year)
 *      Original version                                    01-11-2006
 *      1. Create version for 64bit integers                27-04-2011
 */

/* Table of CRCs of all 8-bit messages. */
static unsigned long crc_table[256];

/* Flag: has the table been computed? Initially false. */
static int crc_table_computed = 0;

/* Make the table for a fast CRC. */
void make_crc32_table(void)
{
    unsigned long c;

    int n, k;
    for (n = 0; n < 256; n++) {
	c = (unsigned long) n;
	for (k = 0; k < 8; k++) {
	    if (c & 1) {
		c = 0xedb88320L ^ (c >> 1);
	    } else {
		c = c >> 1;
	    }
	}
	crc_table[n] = c;
    }
    crc_table_computed = 1;
}

/*
 Update a running crc with the bytes buf[0..len-1] and return
 the updated crc. The crc should be initialized to zero.
*/
unsigned long update_crc32(unsigned long crc,
	                 unsigned char *buf, long len)
{
    unsigned long   c = crc ^ 0xffffffffL;
    long             n;

    if (!crc_table_computed)
	make_crc32_table();
    for (n = 0; n < len; n++) {
	c = crc_table[(c ^ buf[n]) & 0xff] ^ (c >> 8);
    }
    return c ^ 0xffffffffL;
}

/* Return the CRC of the bytes buf[0..len-1]. */
unsigned long compute_crc32(unsigned char *buf, long len)
{
    return update_crc32(0L, buf, len);
}
