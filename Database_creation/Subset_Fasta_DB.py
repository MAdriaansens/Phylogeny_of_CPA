tinput = '/nesi/nobackup/uc04105/new_databases_May/GTDB_226/Bacteria_GTDB226_protein_May92025.faa'

#og script is from BioPython
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

from Bio import SeqIO

record_iter = SeqIO.parse(open(tinput), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, 7070931)):
    filename = "/nesi/nobackup/uc04105/new_databases_May/GTDB_226/DB/Bacteria_GTDB226_protein_May92025_subset%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
~                                                          
