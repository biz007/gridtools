#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
A toolset for quality control, evlation and processing of GRID-seq libarary.
'''

import os, sys, argparse, re
import numpy as np 
import pandas as pd
import pysam 
import h5py
import sqlite3

from concurrent.futures import ProcessPoolExecutor



class GridLinker(object):
    '''
        Class for GRID-seq raw reads mapped to GRID-linker.
    '''

    def __init__(self, bamfile, h5file, minlen, readname):
        self.bamfile = bamfile
        self.h5file = h5file
        self.minlen = minlen
        self.readname = readname
        self._stats = {}
        

    def fqmates(self):
        with pysam.AlignmentFile(self.bamfile, 'rb', require_index=True, threads=os.cpu_count()) as bam: # pylint: disable=maybe-no-member
            self._stats['counts'] = {'reads.total':0, 'reads.linked':0, 'reads.mates':0}
            self._stats['sizes'] = np.zeros((100,3), dtype=int)

            lenr, lend, lenl = 20, 20, bam.lengths[0]
            self._stats['bases'] = np.zeros((lenr+lenl+lend, 5), dtype=int)
            self._stats['quals'] = np.zeros(lenr+lenl+lend, dtype=int)
        
            i = 0
            for read in bam:
                if read.flag == 0 or read.flag == 16:
                    self._stats['counts']['reads.total'] += 1

                    llen = read.reference_length
                    rlen = read.query_alignment_start
                    dlen = read.query_length-read.query_alignment_end
                    self._stats['sizes'][[rlen, llen, dlen], [0,1,2]] += 1

                    if llen == lenl:
                        self._stats['counts']['reads.linked'] += 1

                        if rlen >= self.minlen and dlen >= self.minlen:
                            i += 1
                            self._stats['counts']['reads.mates'] += 1

                            lps = read.query_alignment_start
                            lpe = read.query_alignment_end
                            rid = read.query_name
                            if self.readname:
                                rid = '{}:{}'.format(self.readname, i)

                            rseq = read.query_sequence[0:lps]
                            lseq = read.query_sequence[lps:lpe]
                            dseq = read.query_sequence[lpe:]
                            rqul = read.query_qualities[0:lps]
                            lqul = read.query_qualities[lps:lpe]
                            dqul = read.query_qualities[lpe:]

                            rseq1 = ('N'*lenr + rseq)[-lenr:]
                            dseq1 = (dseq + 'N'*lend)[:lend]
                            lseq1 = (lseq + 'N'*lenl)[:lenl]
                            rqul1 = np.append(np.zeros(lenr, dtype=int), rqul)[-lenr:]
                            dqul1 = np.append(dqul, np.zeros(lend, dtype=int))[:lend]
                            lqul1 = np.append(lqul, np.zeros(lenl, dtype=int))[:lenl]

                            self._stats['quals'] += np.concatenate((rqul1, lqul1, dqul1))
                            
                            base = np.array(list((rseq1 + lseq1 + dseq1).upper()))
                            self._stats['bases'] += np.concatenate(
                                (base == 'A', base == 'C', base == 'G', base == 'T', base == 'N'), 
                                axis = 0).reshape(5, -1).T

                            mstr = '@{}\n{}\n+\n{}\n@{}\n{}\n+\n{}'.format(
                                rid, rseq, ''.join([chr(c+33) for c in rqul]),
                                rid, dseq, ''.join([chr(c+33) for c in dqul])
                            )
                            
                            yield mstr


    def tohdf5(self):
        with h5py.File(self.h5file, 'w', libver='latest') as hf:
            hg = hf.create_group('/linker/stats')
            for k in self._stats.keys():
                if k == 'counts':
                    xa = np.fromiter(self._stats[k].items(), dtype=[('key', 'S20'), ('value', 'i4')])
                    hg.create_dataset(k, data=xa)
                elif k == 'quals':
                    xa = self._stats[k] // self._stats['counts']['reads.mates']
                    hg.create_dataset(k, data=xa)
                else:
                    hg.create_dataset(k, data=self._stats[k], compression='gzip')

    @property
    def stats(self):
        return self._stats

    @property
    def minlen(self):
        return self._minlen

    @minlen.setter
    def minlen(self, minlen):
        if minlen > 0:
            self._minlen = minlen
        else:
            self._minlen = None

    @property
    def readname(self):
        return self._readname

    @readname.setter
    def readname(self, readname):
        if isinstance(readname, str):
            self._readname = readname
        else:
            self._readname = None


    @property
    def bamfile(self):
        return self._bamfile
    
    @bamfile.setter
    def bamfile(self, file):
        self._bamfile = file


    @property
    def h5file(self):
        return self._h5file

    @h5file.setter
    def h5file(self, filename):
        self._h5file = filename


#  ----------------------------------- END of GridLinker ----------------------------------- 


class GridGenome(object):
    '''
        Class for GRID-seq mate reads mapped to the genome.
    '''

    def __init__(self, bamfile, h5file, bink, winm, gtf):
        self.bamfile = bamfile
        self.h5file = h5file
        self.bink = bink
        self.winm = winm
        self.gtf = gtf
        self._stats = {}
        self._data = {}

        ## init files
        self.dfgene
        self.dfchrom


    @property
    def bamfile(self):
        return self._bamfile
    
    @bamfile.setter
    def bamfile(self, file):
        self._bamfile = file

    @property
    def h5file(self):
        return self._h5file

    @h5file.setter
    def h5file(self, filename):
        self._h5file = filename

    @property
    def bink(self):
        return self._bink

    @bink.setter
    def bink(self, k):
        self._bink = k

    @property
    def winm(self):
        return self._winm
        
    @winm.setter
    def winm(self, m):
        self._winm = m

    @property
    def gtf(self):
        return self._gtf

    @gtf.setter
    def gtf(self, fn):
        self._gtf = fn

    @property
    def dfgene(self):
        if 'dfgene' not in self._data:
            xgene = {}
            with pysam.TabixFile(self.gtf) as tbx:                                      # pylint: disable=maybe-no-member
                for gtf in tbx.fetch(parser=pysam.asGTF()):                             # pylint: disable=maybe-no-member
                    if gtf.feature == 'gene':
                        xgene[gtf.gene_id] = (gtf.contig, gtf.start, gtf.end, gtf.strand, gtf.gene_id)
            self._data['dfgene'] = pd.DataFrame(
                np.array(list(xgene.values()), dtype=[
                    ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), 
                    ('Strand', 'S1'), ('GeneID', 'S30')]
                )
            )

        return self._data['dfgene']

    @property
    def dfchrom(self):
        if 'dfchrom' not in self._data:
            with pysam.AlignmentFile(self.bamfile, 'rb') as bam:                        # pylint: disable=maybe-no-member
                self._data['dfchrom'] = pd.DataFrame({
                    'Chrom': [x.encode('ascii') for x in bam.references],
                    'Size': bam.lengths
                })

        return self._data['dfchrom']

    @property
    def dfbinid(self):
        if 'dfbinid' not in self._data:
            chrbins = self.dfchrom['Size'] // (self.bink*1000) + 1
            xs = np.repeat(self.dfchrom['Chrom'], chrbins)
            xr = np.concatenate([np.arange(n) for n in chrbins])
            df = pd.DataFrame({'Chrom': xs, 'Bin': xr})
            df.Bin = df.Bin * self.bink*1000
            
            self._data['dfbinid'] = df
        
        return self._data['dfbinid']


    def itermate(self):
        self._stats['counts'] = {'mates.paired':0, 'mates.mapped':0, 'mates.unique':0, 'mates.dups':0}

        with pysam.AlignmentFile(self.bamfile, 'rb', threads=os.cpu_count()) as bam:                            # pylint: disable=maybe-no-member
            if bam.header['HD']['SO'] != 'queryname':
                raise IOError('The input bam file is not sorted by read name.')

            # rid, RNA [chrom, pos, strand], DNA [chrom, pos, strand], rdup
            mread = {}
            for read in bam:
                if read.is_paired:
                    rid = read.query_name
                    rchr = '.' if read.is_unmapped else read.reference_name
                    rpos = -1 if read.is_unmapped else int((read.reference_start + read.reference_end)/2)
                    rstr = '.' if read.is_unmapped else ('-' if read.is_reverse else '+')
                    rdup = 1 if read.is_duplicate else 0

                    if read.is_read1:
                        mread['rid'] = rid
                        mread['RNA'] = (rchr, rpos, rstr)
                        mread['DNA'] = ('.', -2, '.')
                        mread['rdup'] = rdup
                    elif mread['rid'] == rid:
                        mread['DNA'] = (rchr, rpos, rstr)
                        mread['rdup'] += rdup
                        self._stats['counts']['mates.paired'] += 1

                        if mread['RNA'][1] >= 0 and mread['DNA'][1] >= 0:
                            self._stats['counts']['mates.mapped'] += 1

                            if mread['rdup'] > 1:
                                self._stats['counts']['mates.dups'] += 1
                            else:
                                self._stats['counts']['mates.unique'] += 1
                                yield (mread['rid'],) + mread['RNA'] + mread['DNA']


    def _p_join_dfgene(self, dfx):
        conn = sqlite3.connect(':memory:')
        self.dfgene.to_sql('gene', conn, index=False, index_label=['Chrom', 'Start', 'End', 'Strand'])
        dfx.to_sql('rmate', conn, index=False, index_label=['rchrom', 'rpos', 'rstrand', 'dchrom', 'dpos'])

        qstr = '''
            SELECT 
                seqid, rchrom, rpos, rstrand, dchrom, dpos, dstrand,
                GeneID, CAST(dpos / {0}  as INTEGER) * {0} Bin,
                CASE dchrom = rchrom WHEN 1 THEN 'cis' ELSE 'trans' END xtype
            FROM rmate LEFT JOIN gene ON 
                rchrom = Chrom and rstrand = Strand and rpos BETWEEN Start AND End
        '''.format(self.bink*1000)

        dfy = pd.read_sql_query(qstr, conn)
        conn.close()

        return dfy


    def evalmate(self):
        '''
            store the read mate in dataframe
        '''

        if 'rmate' not in self._data:
            dfrmate = pd.DataFrame(
                np.array(list(self.itermate()), dtype=[
                    ('seqid', 'S100'), 
                    ('rchrom', 'S10'), ('rpos', 'i4'), ('rstrand', 'S1'),
                    ('dchrom', 'S10'), ('dpos', 'i4'), ('dstrand', 'S1')
                ])
            )

            dflist = np.array_split(dfrmate, os.cpu_count()*10, axis=0)

            with ProcessPoolExecutor() as pool:
                self._data['rmate'] = pd.concat(
                    pool.map(self._p_join_dfgene, dflist)
                )

        return True


    def evaldna(self):
        '''
            summarize DNA reads in genomic bins
        '''

        df = self._data['rmate'][
                ['dchrom', 'Bin', 'xtype', 'GeneID']
            ].groupby(['dchrom', 'Bin', 'xtype']).count()['GeneID'].unstack('xtype', fill_value=0).reset_index().rename(
                columns={'dchrom':'Chrom'}
            )
        
        df = pd.merge(self.dfbinid, df, how='left', on=['Chrom', 'Bin']).fillna(0)

        self._data['DNA.reads'] = df

        return True


    def evalrna(self):
        '''
            summarize RNA reads in genes
        '''

        dfscope = pd.merge(
            self.dfgene[['GeneID', 'Chrom', 'Start', 'End', 'Strand']], 
            self._data['rmate'][['GeneID', 'xtype', 'dpos']],
            how='inner', on=['GeneID']
            ).assign(
                Scope = lambda x: np.where(
                    x.xtype == 'trans', -1, 
                    np.int8(np.log10(np.maximum(np.abs(x.dpos - (x.Start+x.End)/2) - (x.End-x.Start)/2, 0)+1))
                )
            ).groupby(['GeneID', 'Chrom', 'Start', 'End', 'Strand', 'Scope']).size().reset_index().rename(columns={0:'reads'})
        
        self._data['RNA.scope'] = dfscope
        
        self._data['RNA.exprs'] = dfscope.groupby(
            ['GeneID', 'Chrom', 'Start', 'End', 'Strand']
        )['reads'].sum().reset_index().assign(
            RPK = lambda x: x.reads/(x.End - x.Start)*1e3,
            TPM = lambda x: x.RPK/np.sum(x.RPK)*1e6
        ).sort_values(by='RPK', ascending=False)

        return True


    def evalmatrix(self):
        '''
            summarize RNA reads in gene-Bin wise
        '''

        self._data['matrix'] = pd.merge(
            self.dfgene[['GeneID', 'Chrom', 'Start', 'End', 'Strand']], 
            self._data['rmate'][['GeneID', 'dchrom', 'Bin']],
            how='inner', on=['GeneID']
            ).groupby(['GeneID', 'Chrom', 'Start', 'End', 'Strand', 'dchrom', 'Bin']).size().reset_index().rename(columns={0:'reads'})

        return True


    def evaluate(self):
        self.evalmate()
        self.evaldna()
        self.evalrna()
        self.evalmatrix()

        return True

    
    def tohdf5(self):
        with h5py.File(self.h5file, 'r+') as hf:
            h5path = '/genome/stats'

            if h5path in hf: del hf[h5path]
            hg = hf.create_group(h5path)
            hg.create_dataset('bink', data=self.bink)
            hg.create_dataset('winm', data=self.winm)
            xa = np.fromiter(self._stats['counts'].items(), dtype=[('key', 'S20'), ('value', 'i4')])
            hg.create_dataset('counts', data=xa)


            h5path = '/genome/files'
            if h5path in hf: del hf[h5path]
            hg = hf.create_group(h5path)
            chrom = self.dfchrom.to_records(index=False).astype(
                [('Chrom', 'S10'), ('Size', 'i4')]
            )
            hg.create_dataset('chrom', data=chrom, compression='gzip')
            
            gene = self.dfgene.to_records(index=False).astype(
                [('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'), ('GeneID', 'S30')]
            )
            hg.create_dataset('gene', data=gene, compression='gzip')


            h5path = '/genome/data'
            if h5path in hf: del hf[h5path]
            hg = hf.create_group(h5path)

            rmate = self._data['rmate'].to_records(index=False).astype(
                [('seqid', 'S100'), 
                ('rchrom', 'S10'), ('rpos', 'i4'), ('rstrand', 'S1'),
                ('dchrom', 'S10'), ('dpos', 'i4'), ('dstrand', 'S1'),
                ('GeneID', 'S30'), ('Bin', 'i4'), ('xtype', 'S6')]
            )
            hg.create_dataset('rmate', data=rmate, compression='gzip')

            matrix = self._data['matrix'].to_records(index=False).astype(
                [('GeneID', 'S30'), ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'),
                ('dchrom', 'S10'), ('Bin', 'i4'), ('reads', 'i4')]
            )
            hg.create_dataset('matrix', data=matrix, compression='gzip')

            gexpr = self._data['RNA.exprs'].to_records(index=False).astype(
                [('GeneID', 'S30'), ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'), 
                ('reads', 'i4'), ('RPK', 'f4'), ('TPM', 'f4')]
            )
            hg.create_dataset('RNA.exprs', data=gexpr, compression='gzip')

            gscop = self._data['RNA.scope'].to_records(index=False).astype(
                [('GeneID', 'S30'), ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'), 
                ('Scope', 'i1'), ('reads', 'i4')]
            )
            hg.create_dataset('RNA.scope', data=gscop, compression='gzip')

            dbins = self._data['DNA.reads'].to_records(index=False).astype(
                [('Chrom', 'S10'), ('Bin', 'i4'), ('cis', 'i4'), ('trans', 'i4')]
            )
            hg.create_dataset('DNA.reads', data=dbins, compression='gzip')

        return True


#  ----------------------------------- END of GridGenome ----------------------------------- 


class GridInfo(object):
    def __init__(self, h5file):
        self.h5file = h5file
        self._bink = None
        self._winm = None
        self._stats = {}
        self._data = {}

    @property
    def h5file(self):
        return self._h5file

    @h5file.setter
    def h5file(self, filename):
        self._h5file = filename

    @property
    def bink(self):
        if not self._bink:
            self._bink = np.int(self._fromhdf5('/data/DNA/bink'))

        return self._bink

    @property
    def winm(self):
        if not self._winm:
            self._winm = np.int(self._fromhdf5('/data/DNA/winm'))

        return self._winm

    @property
    def dfgene(self):
        if 'dfgene' not in self._data:
            self._data['dfgene'] = pd.DataFrame(self._fromhdf5('/files/gene'))

        return self._data['dfgene']

    @property
    def dfchrom(self):
        if 'dfchrom' not in self._data:
            self._data['dfchrom'] = pd.DataFrame(self._fromhdf5('/files/chrom'))

        return self._data['dfchrom']

    def filterChrom(self, pattern='chr([0-9]+|[XY])'):
        ptn = re.compile(pattern)

        chrs = [x for x in self.dfchrom.Chrom if ptn.search(str(x))]
        self._data['dfchrom'] = self.dfchrom[self.dfchrom.Chrom.isin(chrs)]


    def _bpcc(self, data, bink, kfold):
        np.random.shuffle(data)
        xd = [
            np.bincount(x // (bink*1000)) for x in np.array_split(data, kfold)[:kfold]
        ]

        n = min([len(x) for x in xd])
        xd = np.array([x[:n] for x in xd])

        pcc = np.corrcoef(xd)[np.triu_indices(kfold, 1)]
        return np.mean(pcc)


    def _kpcc(self, npa, kbins, kfold=2, nrep=10):
        df = pd.DataFrame(np.array([
            (bink, i, self._bpcc(npa, bink, kfold)) for bink in kbins for i in range(nrep)
        ], dtype=[('BinK', 'i2'), ('Repeat', 'i2'), ('PCC', 'f2')]))

        df = df.set_index('BinK').groupby('BinK')['PCC'].mean()

        return df


    def getResolution(self, kbins):
        df = self.dfrmate[self.dfrmate['rchrom'] == self.dfrmate['dchrom']][['dchrom', 'dpos']]
        df.columns = ['Chrom', 'Pos']
        

        dfpcc = df.set_index('Chrom').groupby('Chrom').apply(
            lambda x: self._kpcc(x.Pos.values, kbins)
        ).reset_index()

        return dfpcc


    @property
    def dfrmate(self):
        if 'rmate' not in self._data:
            df = pd.DataFrame(self._fromhdf5('/data/rmate'))
            df = df[df.rchrom.isin(self.dfchrom.Chrom) & df.dchrom.isin(self.dfchrom.Chrom)]
            self._data['rmate'] = df

        return self._data['rmate']


    @property
    def dfRawMatrix(self):
        if 'rawmtx' not in self._data:
            df = pd.DataFrame(self._fromhdf5('/data/matrix'))
            df = df[df.Chrom.isin(self.dfchrom.Chrom)]
            self._data['rawmtx'] = df

        return self._data['rawmtx']

    @property
    def dfBinReads(self):
        if 'DNA' not in self._data:
            self._data['DNA'] = {}
        if 'reads' not in self._data['DNA']:
            df = pd.DataFrame(self._fromhdf5('/data/DNA/reads'))
            df = df[df.Chrom.isin(self.dfchrom.Chrom)]
            self._data['DNA']['reads'] = df

        return self._data['DNA']['reads']

    @property
    def dfBinValues(self):
        if 'DNA' not in self._data:
            self._data['DNA'] = {}
        if 'values' not in self._data['DNA']:
            df = self.dfBinReads.set_index(['Chrom', 'Bin']).transform(
                lambda x: x/np.sum(x)*1e6
            ).groupby(['Chrom'], sort=False).transform(
                lambda x: np.log2(x.rolling(self.winm*10-1, center=True, min_periods=1).mean()+1)
            ).reset_index()
            self._data['DNA']['values'] = df

        return self._data['DNA']['values']

    def readsToValue(self, df):
        df = df.transform(
            lambda x: np.where(x.rolling(self.winm*2-1, center=True, min_periods=1).sum()>=0.3*self.winm, x, 0)
        ).transform(
            lambda x: x/np.sum(x)*1e6
        ).groupby(['Chrom'], sort=False).transform(
            lambda x: np.log2(x.rolling(self.winm*2-1, center=True, min_periods=1).mean()+1)
        )

        return df

    @property
    def dfGeneExprs(self):
        if 'RNA' not in self._data:
            self._data['RNA'] = {}
        if 'exprs' not in self._data['RNA']:
            df = pd.DataFrame(self._fromhdf5('/data/RNA/exprs'))
            dfd = self.dfRawMatrix.set_index('GeneID').groupby('GeneID')['reads'].max()/self.bink
            df = pd.merge(df, dfd.reset_index().rename(columns={'reads': 'dRPK'}), on='GeneID')
            self._data['RNA']['exprs'] = df

        return self._data['RNA']['exprs']

    @property
    def dfGeneScope(self):
        if 'RNA' not in self._data:
            self._data['RNA'] = {}
        if 'scope' not in self._data['RNA']:
            self._data['RNA']['scope'] = pd.DataFrame(self._fromhdf5('/data/RNA/scope'))

        return self._data['RNA']['scope']


    def v4c(self, geneid):
        df = pd.merge(
            self.dfBinValues[['Chrom', 'Bin']],
            self.dfRawMatrix[self.dfRawMatrix.GeneID == geneid][['Chrom', 'Bin', 'reads']],
            on=['Chrom', 'Bin'], how='left'
        ).set_index(['Chrom', 'Bin']).fillna(0)
        df = self.readsToValue(df).sub(self.dfBinValues.set_index(['Chrom', 'Bin'])['BG'], axis='index').rename(
            columns={'reads': 'V'}
        ).transform(
            lambda x: np.where(x > 0, np.exp2(x), 0)
        ).reset_index().assign(GeneID = geneid)

        return df[['GeneID', 'Chrom', 'Bin', 'V']]

    
    def matrix(self, grpk=100, drpk=10):
        gids = self.dfGeneExprs[(self.dfGeneExprs.RPK >= grpk) & (self.dfGeneExprs.dRPK >= drpk)]['GeneID'].values

        for g in gids:
            yield self.v4c(g)


    def _fromhdf5(self, path):
        with h5py.File(self.h5file, 'r') as hf:
            return np.array(hf.get(path))
                
#  ----------------------------------- GridInfo ----------------------------------- 




#############################################################################################################

if __name__ == '__main__':
    pargs = argparse.ArgumentParser(
        description='Tools for processing and evaluation of GRID-seq library.'
    )
    pargs.add_argument('--version', action='version', version='Version: %(prog)s 1.1')
    subpargs = pargs.add_subparsers(dest='cmd')
    
    pargs_mate = subpargs.add_parser('matefq', help='parse bamfile to RNA/DNA mate fastq.')
    pargs_mate.add_argument('bam', help='BAM file mapped to GRID Linker.')
    pargs_mate.add_argument('-o', '--hdf5', required=True, help='output to hdf5 file.')
    pargs_mate.add_argument('-l', '--minlen', required=False, type = int, default=19, help='mates with both RNA and DNA length >= MINLEN [default:19].')
    pargs_mate.add_argument('-n', '--readname', required=False, help='prefix of mate name; default: the orignial read name.')
    

    pargs_eval = subpargs.add_parser('evaluate', help='evaluate the mates quality and quanitity from the bamfile mapped to the genome.')
    pargs_eval.add_argument('bam', help='BAM file mapped to the genome.')
    pargs_eval.add_argument('-o', '--hdf5', required=True, help='update the hdf5 file.')
    pargs_eval.add_argument('-k', '--bink', default=10, type=int, help='bin size (kb) of the genome.')
    pargs_eval.add_argument('-m', '--winm', default=10, type=int, help='moving bins for smoothing.')
    pargs_eval.add_argument('-g', '--gtf', required=True, help='gene annotation in GTF format.')


    pargs_rna = subpargs.add_parser('RNA', help='identify chromatin-enriched RNAs.')
    pargs_rna.add_argument('hdf5', help='HDF5 file with mapping information.')
    pargs_rna.add_argument('-e', '--exprs', required=False, help='output the gene expression [default: print].')
    pargs_rna.add_argument('-s', '--scope', required=False, help='output the RNA interaction scope [default: None].')


    pargs_dna = subpargs.add_parser('DNA', help='identify RNA-enriched chromatin regions in background (trans) and foreground (cis).')
    pargs_dna.add_argument('hdf5', help='HDF5 file with mapping information.')


    pargs_matrix = subpargs.add_parser('matrix', help='cacluate the RNA-chromatin interation matrix.')
    pargs_matrix.add_argument('hdf5', help='HDF5 file with mapping information.')
    pargs_matrix.add_argument('-k', '--rpk', required=False, type=float, default=100, help='RPK of RNA reads from gene. [default: 100]')
    pargs_matrix.add_argument('-x', '--drpk', required=False, type=float, default=10, help='RPK of DNA reads at the max interacted by gene. [default: 10]')


    pargs_model = subpargs.add_parser('model', help='model the network of enhancer-promoter proximity.')
    pargs_model.add_argument('hdf5', help='HDF5 file with mapping information.')
    pargs_model.add_argument('-k', '--rpk', required=False, type=float, default=100, help='RPK of RNA reads from gene. [default: 100]')
    pargs_model.add_argument('-x', '--drpk', required=False, type=float, default=10, help='RPK of DNA reads at the max interacted by gene. [default: 10]')
    pargs_model.add_argument('-e', '--elebed', required=True, help='regulatory elements in BED format.')
    pargs_model.add_argument('-z', '--zscore', required=False, type=float, default=-10, help='z-score of significant proximity.')
    

    pargs_stats = subpargs.add_parser('stats', help='statistics of GRID-seq data.')
    pargs_stats.add_argument('hdf5', help='HDF5 file with mapping information.')
    pargs_stats.add_argument('-p', '--prefix', required=True, help='prefix for output file names.')
    pargs_stats.add_argument('-c', '--counts', action='store_true', help='if output the summary of mapping information in read counts.')
    pargs_stats.add_argument('-l', '--lengths', action='store_true', help='if output the distribution of sequence length for RNA, DNA and Linker.')
    pargs_stats.add_argument('-b', '--bases', action='store_true', help='if output the summary of base-position information for RNA, DNA and Linker.')
    pargs_stats.add_argument('-q', '--qualities', action='store_true', help='if output the summary of quality-position information for RNA, DNA and Linker.')
    pargs_stats.add_argument('-r', '--resolution', action='store_true', help='if output the resolution information of the library.')


    args = pargs.parse_args()

    if args.cmd == 'matefq':
        grid = GridLinker(args.bam, args.hdf5, args.minlen, args.readname)

        for fq in grid.fqmates():
            print(fq)

        grid.tohdf5()


    elif args.cmd == 'evaluate':
        grid = GridGenome(args.bam, args.hdf5, args.bink, args.winm, args.gtf)
        grid.evaluate()
        grid.tohdf5()


    elif args.cmd == 'RNA':
        grid = GridInfo(args.hdf5)
        grid.filterChrom()

        dfexpr = grid.dfGeneExprs[['GeneID', 'reads', 'TPM', 'RPK', 'dRPK']]
        dfexpr.GeneID = dfexpr.GeneID.values.astype('U30')

        if args.exprs:
            dfexpr.to_csv(args.exprs, header=True, index=False, sep='\t', float_format='%.3f')
        else:
            print(dfexpr.to_csv(header=True, index=False, sep='\t', float_format='%.3f'))
        
        if args.scope:
            dfscope = grid.dfGeneScope[['GeneID', 'Scope', 'reads']]
            dfscope.GeneID = dfscope.GeneID.values.astype('U30')
            dfscope.to_csv(args.scope, header=True, index=False, sep='\t', float_format='%.3f')


    elif args.cmd == 'DNA':
        grid = GridInfo(args.hdf5)
        grid.filterChrom()

        df = grid.dfBinValues
        df.Chrom = df.Chrom.values.astype('U10')
   
        print(df.to_csv(header=True, index=False, sep='\t', float_format='%.3f'))


    elif args.cmd == 'matrix':
        grid = GridInfo(args.hdf5)
        grid.filterChrom()

        for df in grid.matrix(args.rpk, args.drpk):
            df.GeneID = df.GeneID.values.astype('U30')
            df.Chrom = df.Chrom.values.astype('U10')
            print(df.to_csv(header=True, index=False, sep='\t', float_format='%.3f'))


    elif args.cmd == 'model':
        grid = GridInfo(args.hdf5)
        grid.filterChrom()
        binsize = grid.bink*1e3
        
        conn = sqlite3.connect(':memory:')

        grid.dfgene.to_sql('gene', conn, index=False, index_label=['GeneID', 'Chrom'])

        dfbed = pd.read_csv(args.elebed, sep='\t', header=None, 
            usecols=[0,1,2,3], names=['Chrom', 'Start', 'End', 'EType']).assign(
                EID = lambda x: x.Chrom.str.cat([map(str, x.Start), map(str, x.End)], sep=':')
            )
        dfbed.Chrom = dfbed.Chrom.values.astype('S10')
        dfbed.Start = np.int32(np.floor(dfbed.Start.values / binsize))
        dfbed.End = np.int32(np.ceil(dfbed.End.values / binsize))

        
        dfbed.to_sql('bed', conn, index=False, index_label=['Chrom', 'Start', 'End'])

        qstr = '''
            SELECT 
                U.GeneID GeneID, EID, EType,
                CASE U.Chrom = gene.Chrom WHEN 1 THEN 'cis' ELSE 'trans' END XType,
                S, V
            FROM (
                SELECT 
                    GeneID, mtx.Chrom Chrom, EID, EType, S, SUM(V) V
                FROM
                    mtx JOIN bed ON mtx.Chrom = bed.Chrom AND 
                    mtx.Bin BETWEEN bed.Start AND bed.END
                GROUP BY
                    GeneID, mtx.Chrom, EID, EType, S
            ) U JOIN gene ON U.GeneID = gene.GeneID 
        '''

        for dfmat in grid.matrix(args.rpk, args.drpk):
            dfmat = dfmat[dfmat.V>0]
            dfmat['S'] = np.sum(dfmat.V.values)
            dfmat.to_sql('mtx', conn, index=False, index_label=['GeneID', 'Chrom', 'Bin'], if_exists='append')

        df = pd.read_sql_query(qstr, conn)
        conn.close()
        
        df = df.assign(
            G = lambda x: np.log10(x.V)-np.log10(x.S)
        ).assign(
            Z = lambda x: (x.G - x.G[x.XType == 'trans'].mean())/x.G[x.XType == 'trans'].std()
        )

        df = df[df.Z >= args.zscore]
        df.GeneID = df.GeneID.values.astype('U30')

        print(df.to_csv(header=True, index=False, sep='\t', float_format='%.3f', 
            columns=['GeneID', 'EID', 'EType', 'XType', 'G', 'Z']))


    elif args.cmd == 'stats':
        grid = GridInfo(args.hdf5)
        grid.filterChrom()

        with h5py.File(args.hdf5, 'r') as hf:
            if args.bases:
                df = pd.DataFrame(np.array(hf.get('stats/bases')))
                df.columns = list('ACGTN')
                df['base'] = df.index+1
                df = df[['base', 'A', 'C', 'G', 'T', 'N']]
                df.to_csv(args.prefix + '.bases.txt', header=True, index=False, sep='\t', float_format='%i')
            if args.counts:
                df = pd.DataFrame(np.array(hf.get('stats/counts')))
                df.key = df.key.values.astype('U20')
                df.to_csv(args.prefix + '.counts.txt', header=False, index=False, sep='\t', float_format='%i')
            if args.lengths:
                df = pd.DataFrame(np.array(hf.get('stats/sizes')))
                df.columns = ['RNA', 'linker', 'DNA']
                df['length'] = df.index
                df = df[['length', 'RNA', 'linker', 'DNA']]
                df.to_csv(args.prefix + '.lengths.txt', header=True, index=False, sep='\t', float_format='%i')
            if args.qualities:
                quals = np.array(hf.get('stats/quals'))
                df = pd.DataFrame({
                    'base': np.arange(1, len(quals)+1),
                    'qual': quals })
                df.to_csv(args.prefix + '.quals.txt', header=True, index=False, sep='\t', float_format='%i')
            if args.resolution:
                # kbins = [1,10,100,1000,10000]
                kbins = [1, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 1000]
                df = grid.getResolution(kbins)
                df.Chrom = df.Chrom.values.astype('U10')
                df.columns = ['Chrom']+[str(k)+'K' for k in kbins]
                df.to_csv(args.prefix + '.resolution.txt', header=True, index=False, sep='\t', float_format='%.3f')


    else:
        pargs.print_help()
