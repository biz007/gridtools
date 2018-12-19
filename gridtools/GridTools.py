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
from itertools import islice, chain, starmap



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
        with h5py.File(self.h5file, libver='latest') as hf:
            h5path = '/linker/stats'

            if h5path in hf: del hf[h5path]
            hg = hf.create_group(h5path)
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
        self.gtffile = gtf
        self._stats = {}
        self._data = {}

        ## init files
        self.dfgene
        self.dfchrom

        # with h5py.File(self.h5file, libver='latest') as hf:
        #     h5path = '/genome'
        #     if h5path in hf: del hf[h5path]
        #     hfg = hf.create_group(h5path)
        #     hfg.create_dataset('stats/bink', data=self.bink)
        #     hfg.create_dataset('stats/winm', data=self.winm)

        #     with pysam.AlignmentFile(self.bamfile, 'rb') as bam:                            # pylint: disable=maybe-no-member
        #         ds = pd.DataFrame({
        #             'Chrom': [x.encode('ascii') for x in bam.references],
        #             'Size': bam.lengths
        #         }).to_records(index=False).astype(
        #             [('Chrom', 'S10'), ('Size', 'i4')]
        #         )
                
        #         hfg.create_dataset('files/chrom', data=ds, compression='gzip')

            # with pysam.TabixFile(self.gtffile) as tbx:                                      # pylint: disable=maybe-no-member
            #     glist = []
            #     for gtf in tbx.fetch(parser=pysam.asGTF()):                                 # pylint: disable=maybe-no-member
            #         if gtf.feature == 'gene':
            #             grow =  (gtf.contig, gtf.start, gtf.end, gtf.strand, gtf.gene_id)
            #             glist.append(grow)
                
            #     ds = np.array(glist, dtype=[
            #             ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), 
            #             ('Strand', 'S1'), ('GeneID', 'S30')
            #             ]
            #         )

            #     hfg.create_dataset('files/gene', data=ds, compression='gzip')

            # hf.swmr_mode = True
            # assert hf.swmr_mode


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
    def gtffile(self):
        return self._gtffile

    @gtffile.setter
    def gtffile(self, fn):
        self._gtffile = fn


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
    def dfgene(self):
        if 'dfgene' not in self._data:
            xgene = {}
            with pysam.TabixFile(self.gtffile) as tbx:                                      # pylint: disable=maybe-no-member
                for gtf in tbx.fetch(parser=pysam.asGTF()):                                 # pylint: disable=maybe-no-member
                    if gtf.feature == 'gene':
                        xgene[gtf.gene_id] = (gtf.contig, gtf.start, gtf.end, gtf.strand, gtf.gene_id)
            self._data['dfgene'] = pd.DataFrame(
                np.array(list(xgene.values()), dtype=[
                    ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), 
                    ('Strand', 'S1'), ('GeneID', 'S30')]
                )
            )

        return self._data['dfgene']


    def itermate(self):
        self._stats['counts'] = {'mates.paired':0, 'mates.mapped':0, 'mates.unique':0, 'mates.dups':0}

        with pysam.AlignmentFile(self.bamfile, 'rb') as bam:                            # pylint: disable=maybe-no-member
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


    def iterchunk(self, x, n):
        xlist = list(islice(x, n))

        while xlist:
            yield xlist
            xlist = list(islice(x, n))


    def _group_by_bins(self, dfx):
        df = dfx.assign(
            Bin = lambda df: df.dpos // (self.bink*1000) * 1000,
            xtype = lambda x: np.where(x.dchrom == x.rchrom, b'cis', b'trans')
        ).groupby(['dchrom', 'Bin', 'xtype']).size().reset_index().rename(
            columns={0:'reads'}
        )

        return df


    def _join_by_range(self, dfx):
        with sqlite3.connect(':memory:') as conn:
            self.dfgene.to_sql('gene', conn, index=False, index_label=['Chrom', 'Start', 'End', 'Strand'])
            dfx.to_sql('rmate', conn, index=False, index_label=['rchrom', 'rpos', 'rstrand', 'dchrom', 'dpos'])

            qstr = '''
                SELECT 
                    seqid, rchrom, rpos, rstrand, dchrom, dpos, dstrand, GeneID
                FROM rmate LEFT JOIN gene ON 
                    rchrom = Chrom and rstrand = Strand and rpos BETWEEN Start AND End
            '''

            dfy = pd.read_sql_query(qstr, conn)

        return dfy

    
    def _group_by_gscope(self, dfx):
        dfy = pd.merge(
            self.dfgene[['GeneID', 'Chrom', 'Start', 'End', 'Strand']],
            dfx[['GeneID', 'dchrom', 'dpos']], how='inner', on=['GeneID']
        ).assign(
            Scope = lambda x: np.where(x.Chrom != x.dchrom, -1, 
            np.int8(np.log10(np.maximum(np.abs(x.dpos - (x.Start+x.End)/2) - (x.End-x.Start)/2, 0)+1))
        )).groupby(['GeneID', 'Chrom', 'Start', 'End', 'Strand', 'Scope']).size().reset_index().rename(columns={0:'reads'})

        return dfy
 

    def _group_by_genebins(self, dfx):
        dfy = pd.merge(
            self.dfgene[['GeneID', 'Chrom', 'Start', 'End', 'Strand']],
            dfx.assign(
                Bin = lambda df: df.dpos // (self.bink*1000) * 1000
            ).groupby(['GeneID', 'dchrom', 'Bin']).size().reset_index().rename(columns={0:'reads'})
        )

        return dfy


    def evaluate(self):
        '''
            store the read mate in dataframe
        '''

        dfdna = []; dfmate = []; dfrscope = []; dfgbmtx = []
        for mlist in self.iterchunk(self.itermate(), 1000000):
            df = pd.DataFrame(
                np.array(list(mlist), dtype=[
                ('seqid', 'S50'), 
                ('rchrom', 'S10'), ('rpos', 'i4'), ('rstrand', 'S1'),
                ('dchrom', 'S10'), ('dpos', 'i4'), ('dstrand', 'S1')])
            )

            dflist = np.array_split(df, os.cpu_count()*3, axis=0)
            with ProcessPoolExecutor() as pool:
                dfdna.append(pool.map(self._group_by_bins, dflist))
                dfm = pd.concat(pool.map(self._join_by_range, dflist))
                dfmate.append(dfm)

                dmlist = np.array_split(dfm, os.cpu_count()*3, axis=0)
                dfrscope.append(pool.map(self._group_by_gscope, dmlist))
                dfgbmtx.append(pool.map(self._group_by_genebins, dmlist))

        dfmate = pd.concat(dfmate)
        dfdna = pd.concat(chain(*dfdna)).groupby(['dchrom', 'Bin', 'xtype'])['reads'].sum().reset_index()
        dfrscope = pd.concat(chain(*dfrscope)).groupby(['GeneID', 'Chrom', 'Start', 'End', 'Strand', 'Scope'])['reads'].sum().reset_index()
        dfgbmtx = pd.concat(chain(*dfgbmtx)).groupby(['GeneID', 'Chrom', 'Start', 'End', 'Strand', 'dchrom', 'Bin'])['reads'].sum().reset_index()

        dfrexprs = dfgbmtx.groupby(['GeneID', 'Chrom', 'Start', 'End', 'Strand'])['reads'].agg(['sum', 'max']).rename(
            columns={'sum':'reads', 'max':'dreads'}
        ).reset_index().assign(
            RPK = lambda x: x.reads/(x.End - x.Start)*1e3,
            TPM = lambda x: x.RPK/np.sum(x.RPK)*1e6,
            dRPK = lambda x: x.dreads/self.bink
        ).sort_values(by='RPK', ascending=False).drop(columns=['dreads'])
        

        with h5py.File(self.h5file, libver='latest') as hf:
            h5path = '/genome'
            if h5path in hf: del hf[h5path]
            hfg = hf.create_group(h5path)

            hfg.create_dataset('data/DNA.reads', compression='lzf', 
                data = dfdna.to_records(index=False).astype(
                    [('Chrom', 'S10'), ('Bin', 'i4'), ('xtype', 'S5'), ('reads', 'i4')]
                )
            )

            hfg.create_dataset('data/RNA.scope', compression='lzf', 
                data = dfrscope.to_records(index=False).astype(
                    [('GeneID', 'S30'), ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'), 
                    ('Scope', 'i1'), ('reads', 'i4')]
                )
            )

            hfg.create_dataset('data/RNA.exprs', compression='lzf', 
                data = dfrexprs.to_records(index=False).astype(
                    [('GeneID', 'S30'), ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'), 
                    ('reads', 'i4'), ('RPK', 'f4'), ('TPM', 'f4'), ('dRPK', 'f4')]
                )
            )

            hfg.create_dataset('data/rmate', compression='lzf', 
                data = dfmate.to_records(index=False).astype(
                    [
                        ('seqid', 'S100'), 
                        ('rchrom', 'S10'), ('rpos', 'i4'), ('rstrand', 'S1'),
                        ('dchrom', 'S10'), ('dpos', 'i4'), ('dstrand', 'S1'),
                        ('GeneID', 'S30')
                    ]
                )
            )

            hfg.create_dataset('data/matrix', compression='lzf', 
                data = dfgbmtx.to_records(index=False).astype(
                    [
                        ('GeneID', 'S30'), ('Chrom', 'S10'), ('Start', 'i4'), ('End', 'i4'), ('Strand', 'S1'),
                        ('dchrom', 'S10'), ('Bin', 'i4'), ('reads', 'i4')
                    ]
                )
            )

            xa = np.fromiter(self._stats['counts'].items(), dtype=[('key', 'S20'), ('value', 'i4')])
            hfg.create_dataset('stats/counts', data=xa)
            hfg.create_dataset('stats/bink', data=self.bink)
            hfg.create_dataset('stats/winm', data=self.winm)

            hfg.create_dataset('files/chrom', compression='lzf', 
                data = self.dfchrom.to_records(index=False).astype(
                    [('Chrom', 'S10'), ('Size', 'i4')]
                )
            )
                    
        return True


#  ----------------------------------- END of GridGenome ----------------------------------- 


class GridInfo(object):
    def __init__(self, h5file):
        self.h5file = h5file
        self._data = {}

        with h5py.File(self.h5file, 'r') as hf:
            self._bink = np.array(hf['/genome/stats/bink'])
            self._winm = np.array(hf['/genome/stats/winm'])
            self._dfchrom = pd.DataFrame(
                np.array(hf['/genome/files/chrom'])
            )

        chrs = [x for x in self._dfchrom.Chrom.values if re.search('[cC]hr[^CM].*', str(x))]
        self._dfchrom = self._dfchrom.assign(
            main = lambda x: x.Chrom.isin(chrs)
        )


    @property
    def h5file(self):
        return self._h5file

    @h5file.setter
    def h5file(self, filename):
        self._h5file = filename

    @property
    def bink(self):
        return self._bink

    @property
    def winm(self):
        return self._winm

    @property
    def dfchrom(self):
        return self._dfchrom


    def _bpcc(self, data, bink, kfold):
        np.random.shuffle(data)
        xd = [
            np.bincount(np.maximum(1, x // (bink*1000))) for x in np.array_split(data, kfold)[:kfold]
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
        with h5py.File(self.h5file, 'r') as hf:
            df = pd.DataFrame(np.array(hf['/genome/data/RNA.exprs']))

        df = df[df.rchrom == df.dchrom][['dchrom', 'dpos']].rename(
            columns={'dchrom':'Chrom', 'dpos':'Pos'}
        )
        
        dfpcc = df.set_index('Chrom').groupby('Chrom').apply(
            lambda x: self._kpcc(x.Pos.values, kbins)
        ).reset_index()

        return dfpcc


    @property
    def dfBinReads(self):
        if 'DNA.reads' not in self._data:
            df = self.dfchrom[self.dfchrom.main].assign(
                nBin = lambda x: x.Size.values // (self.bink*1000) + 1
            )[['Chrom', 'nBin']]

            df = pd.DataFrame({
                'Chrom': np.repeat(df.Chrom.values, df.nBin.values),
                'Bin': np.concatenate([np.arange(n)*1000 for n in df.nBin.values])
            })

            with h5py.File(self.h5file, 'r') as hf:
                dfdna = pd.DataFrame(np.array(hf['/genome/data/DNA.reads']))

            df = pd.merge(
                df, dfdna.set_index(['Chrom', 'Bin', 'xtype'])['reads'].unstack('xtype', fill_value=0).reset_index(),
                how='left', on=['Chrom', 'Bin']
            ).rename(columns={b'cis':'FG', b'trans':'BG'})

            self._data['DNA.reads'] = df

        return self._data['DNA.reads']


    @property
    def dfBinValues(self):
        if 'DNA.values' not in self._data:
            df = self.dfBinReads.set_index(['Chrom', 'Bin']).transform(
                lambda x: x/np.sum(x)*1e6
            ).groupby(['Chrom'], sort=False).transform(
                #TODO: evaluate performance for: 
                #lambda x: (x.ewm(span=5).mean() + x[::-1].ewm(span=5).mean()[::-1])/2 
                lambda x: np.log2(x.rolling(self.winm*10-1, center=True, min_periods=1).mean()+1)
            ).reset_index()
            self._data['DNA.values'] = df

        return self._data['DNA.values']


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
        if 'RNA.exprs' not in self._data:
            with h5py.File(self.h5file, 'r') as hf:
                self._data['RNA.exprs']=pd.DataFrame(np.array(hf['/genome/data/RNA.exprs']))

        return self._data['RNA.exprs']


    @property
    def dfGeneScope(self):
        if 'RNA.scope' not in self._data:
            with h5py.File(self.h5file, 'r') as hf:
                self._data['RNA.scope']=pd.DataFrame(np.array(hf['/genome/data/RNA.scope']))
        
        return self._data['RNA.scope']


    @property
    def dfMatrix(self):
        if 'rawmtx' not in self._data:
            with h5py.File(self.h5file, 'r') as hf:
                df=pd.DataFrame(np.array(hf['/genome/data/matrix']))
            
            mainchroms = self.dfchrom.Chrom[self.dfchrom.main].values
            self._data['rawmtx'] = df[df.dchrom.isin(mainchroms)]

        return self._data['rawmtx']


    def v4c(self, geneid, zerobin=True):
        df = pd.merge(
            self.dfBinValues[['Chrom', 'Bin']],
            self.dfMatrix[self.dfMatrix.GeneID == geneid][['dchrom', 'Bin', 'reads']].rename(columns={'dchrom':'Chrom'}),
            on=['Chrom', 'Bin'], how='left'
        ).set_index(['Chrom', 'Bin']).fillna(0)

        dfbg = self.dfBinValues.set_index(['Chrom', 'Bin'])['BG']

        df = self.readsToValue(df)['reads'].sub(dfbg, axis='index').transform(
            lambda x: np.where(x > 0, np.exp2(x), 0)
        ).reset_index().rename(columns={0: 'V'}).assign(GeneID = geneid)

        df = df[['GeneID', 'Chrom', 'Bin', 'V']]
        if not zerobin:
            df = df[df.V > 0]

        return df


    def hitGeneID(self, grpk=100, drpk=10):
        gids = self.dfGeneExprs[
            (self.dfGeneExprs.RPK >= grpk) & (self.dfGeneExprs.dRPK >= drpk)
            ]['GeneID'].values

        return gids

    
    def iterMatrix(self, geneids, zerobin=True, chunks=1):
        mlist = []; i = 0
        for g in geneids:
            i += 1
            if i % chunks:
                mlist.append(self.v4c(g, zerobin))
            else:
                yield mlist
                mlist = []
        yield mlist


    def model(self, bedfile, grpk=100, drpk=10):
        dfbed = pd.read_csv(bedfile, sep='\t', header=None, 
            usecols=[0,1,2,3], names=['Chrom', 'Start', 'End', 'EType']).assign(
                EID = lambda x: x.Chrom.str.cat([map(str, x.Start), map(str, x.End)], sep=':')
            )
        dfbed.Chrom = dfbed.Chrom.values.astype('S10')
        dfbed.Start = np.int32(np.floor(dfbed.Start.values / (grid.bink*1000))) * 1000
        dfbed.End = np.int32(np.ceil(dfbed.End.values / (grid.bink*1000))) * 1000

        for xlist in self.iterchunk(self.matrix(grpk, drpk), 1000):

            with ProcessPoolExecutor() as pool:
                dfdna.append(pool.map(self._group_by_bins, xlist))

            #TODO





        return None
#  ----------------------------------- GridInfo ----------------------------------- 




#############################################################################################################

if __name__ == '__main__':
    pargs = argparse.ArgumentParser(
        description='Tools for processing and evaluation of GRID-seq library.'
    )
    pargs.add_argument('--version', action='version', version='Version: %(prog)s 1.2')
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


    elif args.cmd == 'RNA':
        grid = GridInfo(args.hdf5)

        if args.exprs:
            dfexpr = grid.dfGeneExprs[['GeneID', 'reads', 'TPM', 'RPK', 'dRPK']].assign(
                GeneID = lambda df: df.GeneID.values.astype('U30')
            )
            dfexpr.to_csv(args.exprs, header=True, index=False, sep='\t', float_format='%.3f')
        
        if args.scope:
            dfscope = grid.dfGeneScope[['GeneID', 'Scope', 'reads']].assign(
                GeneID = lambda df: df.GeneID.values.astype('U30')
            )
            dfscope.to_csv(args.scope, header=True, index=False, sep='\t', float_format='%.3f')


    elif args.cmd == 'DNA':
        grid = GridInfo(args.hdf5)

        df = grid.dfBinValues.assign(
            Chrom = lambda df: df.Chrom.values.astype('U10')
        )
   
        print(df.to_csv(header=True, index=False, sep='\t', float_format='%.3f'))


    elif args.cmd == 'matrix':
        grid = GridInfo(args.hdf5)

        gids = grid.hitGeneID(args.rpk, args.drpk)

        for df in grid.iterMatrix(gids, chunks=10):
            df = pd.concat(df).assign(
                GeneID = lambda df: df.GeneID.values.astype('U30'),
                Chrom = lambda df: df.Chrom.values.astype('U10')
            )

            print(df.to_csv(header=True, index=False, sep='\t', float_format='%.3f'))


    elif args.cmd == 'model':
        grid = GridInfo(args.hdf5)

        # conn = sqlite3.connect(':memory:')

        # dfbed = pd.read_csv(args.elebed, sep='\t', header=None, 
        #     usecols=[0,1,2,3], names=['Chrom', 'Start', 'End', 'EType']).assign(
        #         EID = lambda x: x.Chrom.str.cat([map(str, x.Start), map(str, x.End)], sep=':')
        #     )
        # dfbed.Chrom = dfbed.Chrom.values.astype('S10')
        # dfbed.Start = np.int32(np.floor(dfbed.Start.values / (grid.bink*1000))) * 1000
        # dfbed.End = np.int32(np.ceil(dfbed.End.values / (grid.bink*1000))) * 1000
        
        # dfbed.to_sql('bed', conn, index=False, index_label=['Chrom', 'Start', 'End'])

        # qstr = '''
        #     SELECT 
        #         U.GeneID GeneID, EID, EType,
        #         CASE U.Chrom = gene.Chrom WHEN 1 THEN 'cis' ELSE 'trans' END XType,
        #         S, V
        #     FROM (
        #         SELECT 
        #             GeneID, mtx.Chrom Chrom, EID, EType, S, SUM(V) V
        #         FROM
        #             mtx JOIN bed ON mtx.Chrom = bed.Chrom AND 
        #             mtx.Bin BETWEEN bed.Start AND bed.END
        #         GROUP BY
        #             GeneID, mtx.Chrom, EID, EType, S
        #     ) U JOIN gene ON U.GeneID = gene.GeneID 
        # '''

        # for dfmat in grid.matrix(args.rpk, args.drpk):
        #     dfmat = dfmat[dfmat.V>0]
        #     dfmat['S'] = np.sum(dfmat.V.values)
        #     dfmat.to_sql('mtx', conn, index=False, index_label=['GeneID', 'Chrom', 'Bin'], if_exists='append')

        # df = pd.read_sql_query(qstr, conn)
        # conn.close()
        
        # df = df.assign(
        #     G = lambda x: np.log10(x.V)-np.log10(x.S)
        # ).assign(
        #     Z = lambda x: (x.G - x.G[x.XType == 'trans'].mean())/x.G[x.XType == 'trans'].std()
        # )

        # df = df[df.Z >= args.zscore]
        # df.GeneID = df.GeneID.values.astype('U30')

        # print(df.to_csv(header=True, index=False, sep='\t', float_format='%.3f', 
        #     columns=['GeneID', 'EID', 'EType', 'XType', 'G', 'Z']))


    elif args.cmd == 'stats':
        grid = GridInfo(args.hdf5)

        with h5py.File(args.hdf5, 'r') as hf:
            if args.counts:
                dfl = pd.DataFrame(np.array(hf.get('/linker/stats/counts')))
                dfg = pd.DataFrame(np.array(hf.get('/genome/stats/counts')))

                df = pd.concat([dfl, dfg]).assign(
                    key = lambda df: df.key.values.astype('U20')
                )
                df.to_csv(args.prefix + '.counts.txt', header=False, index=False, sep='\t', float_format='%i')

            if args.bases:
                df = pd.DataFrame(np.array(hf.get('/linker/stats/bases')))
                df.columns = list('ACGTN')
                df['qual'] = np.array(hf.get('/linker/stats/quals'))
                df['base'] = df.index+1
                df = df[['base', 'A', 'C', 'G', 'T', 'N', 'qual']]
                df.to_csv(args.prefix + '.bases.txt', header=True, index=False, sep='\t', float_format='%i')

            if args.lengths:
                df = pd.DataFrame(np.array(hf.get('/linker/stats/sizes')))
                df.columns = ['RNA', 'linker', 'DNA']
                df['length'] = df.index
                df = df[['length', 'RNA', 'linker', 'DNA']]
                df.to_csv(args.prefix + '.lengths.txt', header=True, index=False, sep='\t', float_format='%i')

            if args.resolution:
                # kbins = [1,10,100,1000,10000]
                kbins = [1, 2, 4, 6, 8, 10, 20, 40, 60, 80, 100, 1000]
                df = grid.getResolution(kbins)
                df.Chrom = df.Chrom.values.astype('U10')
                df.columns = ['Chrom']+[str(k)+'K' for k in kbins]
                df.to_csv(args.prefix + '.resolution.txt', header=True, index=False, sep='\t', float_format='%.3f')


    else:
        pargs.print_help()
