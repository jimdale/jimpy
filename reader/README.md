Installation Instructions for Jim's sim dump reader
---------------------------------------------------

Use f2py to copmile the readers into python shared object files:

    f2py -c -m fort_dump_read fort_dump_read.f
    f2py -c -m fort_dump_read_endian fort_dump_read_endian.f
