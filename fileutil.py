import gzip, sys, binascii


def is_gz(filename):
    "Tells if the given filename leads to invalid file (-1), gzipped file (0), not non-gzipped file(1)."
    try:
        file = open(filename, 'rb')
    except:
        print('[warning::is_gz] File ({0}) not found. Continue anyway.'.format(filename), file=sys.stderr)
        return -1
    if binascii.hexlify(file.read(2)) == b'1f8b':
        file.close()
        return 0
    else:
        file.close()
        return 1


def opener(filename):
    "Determines if the file is gzipped and deal with it."
    status = is_gz(filename)
    if status==-1:
        print('[error::myopen] File ({0}) not found. Abort.'.format(filename), file=sys.stderr)
        exit(1)
    elif status==0: gzipped = True
    elif status==1: gzipped = False
    if gzipped:
        with gzip.open(filename, 'rb') as file:
            for line in file:
                yield line.decode()
    else:
        with open(filename) as file:
            for line in file:
                yield line
