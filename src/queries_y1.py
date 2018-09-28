"""
This module uses que lib SqlAlchemy to write queries for generate VACs -Value added catalogs-.

An operation represents a query that is built based on the input configuration
-input_params- and optionally, it can depend on intermediate tables -intermediate_tables-.
This way, a single operation can be used to compose many queries.

Intermediate table: is a table created on the DB to store temporary data that 
are used to compute the final result. These tables can either be 'permanent' 
or 'temporary' depending on the kind of operation. 

"""

from sqlalchemy.sql import and_, or_
from sqlalchemy import Table, cast, Integer, func, case
from sqlalchemy.sql.expression import literal_column, between, select

from src import sqlalchemy_extension as se
from db import DbConnection


class Operation:
    """
    Base class to create new operations.

    Basically, when a new operation must be written, we heritage the
Operation class and override the methods validate and select.
    """

    def __init__(self, input_params, intermediate_tables):
        """
            Args:
                input_params: a dictionary that has specific information about the
                operation.
        
                intermediate_tables: It has a list of intermediate tables in which this new
                operation depends.
        """
        self.dal = DbConnection()
        if not self.dal.is_connection_initialized():
            raise("The connection to the database was not initialized")

        self.input_params = input_params
        self.intermediate_tables = intermediate_tables

        self.validate()
        self.stm = self.select()

    def validate(self):
        raise NotImplementedError("Implement this method")

    def select(self):
        """
        This method defines the operation. 

        Returns:  
            select: It must return a SQLAlchemy select statement.
        """
        raise NotImplementedError("Implement this method")

    def sql(self):
        return str(self.stm.compile(compile_kwargs={"literal_binds": True}))

    def create_table(self, table_name):
        self.dal.create_table(table_name, self.stm)

    def delete_table(self, table_name):
        self.dal.create_table(table_name)


class GreatEqualThanSignalColumn(Operation):
    """ Select data where the column signal has value greater or equal to 
    input_params['signal']
    """
    def __init__(self, input_params, intermediate_tables):
        """
            Args:
                dal: database connection class. 

                input_params: a dict containing the keys:
                'table' - specifies the table name to load data
                
                'schema' - specifies the schema name where the table is located
                
                'signal' - value to apply the filter condition
        """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['schema', 'table', 'signal']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input params")

    def select(self):
        table = Table(self.input_params['table'], self.dal.meta, autoload=True,
                      schema=self.input_params['schema'])
        stm = select(
          [table]).where(table.c.signal >= self.input_params['signal'])
        return stm


class CombinedMaps(Operation):
    """ Inner join between a list of tables made on the pixel column.
    """
    def __init__(self, input_params, intermediate_tables):
        """
            Args:
                dal: database connection class. 

                intermediate_tables: a dict containing the keys:
                'tables' - list of dicts in which each dict must have the
                 keys 'table' and 'schema'
        """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['tables']
        if not set(keys).issubset(self.intermediate_tables.keys()):
            raise ValueError("Missing input params")

        if len(self.intermediate_tables['tables']) < 2:
            raise ValueError("At least 2 tables must be declared")

        for table in self.intermediate_tables['tables']:
            if 'table' not in table or 'schema' not in table:
                raise ValueError("each dict must have the keys 'table' and 'schema'")

    def select(self):
        # load tables.
        sub_tables = []
        for table in self.intermediate_tables['tables']:
            sub_tables.append(Table(table['table'], self.dal.meta,
                                    schema=table['schema'], autoload=True))

        # join statement
        stm = select([literal_column("1").label('signal'),
                      sub_tables[0].c.pixel,
                      sub_tables[0].c.ra,
                      sub_tables[0].c.dec])
        stm_join = sub_tables[0]
        for table in sub_tables[1:]:
            stm_join = stm_join.join(table, sub_tables[0].c.pixel ==
                                     table.c.pixel)
        stm = stm.select_from(stm_join)
        return stm


class BadRegions(Operation):
    """ The Bad Regions mask flag regions with artifacts in the release
    """
    def __init__(self, input_params, intermediate_tables):
        """
            Args:
                dal: database connection class. 

                input_params: a dict containing the keys:
                'table' - specifies the table name to load data
                
                'schema' - specifies the schema name where the table is located
                
                'filters' - list of FLAGS to mask regions with artifacts in
                 the release. The meaning fo each FLAG is:
                    1 - Regions with bad astrometric colors
                    2 - Fainter 2MASS star region (8 < J < 12)
                    4 - Large nearby object (R3C catalog)
                    8 - Bright 2MASS star region (5 < J < 8)
                    16 - Near the LMC
                    32 - Yale Bright Star region
                    64 - High density of crazy colors
                    128 - Globular Clusters (William et al. 2010)
        """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['schema', 'table', 'filters']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input params")

        if len(self.input_params['filters']) <= 0:
            raise ValueError("At least one filter condition must be added")

    def select(self):
        table = Table(self.input_params['table'], self.dal.meta, autoload=True,
                      schema=self.input_params['schema'])
        stm = select([table]).where(se.BitwiseAnd(
                                    cast(table.c.signal, Integer),
                                    sum(self.input_params['filters'])) >
                                    literal_column('0'))
        return stm


class Footprint(Operation):
    """ The Footprint combines the Detection Fraction maps, the Bad
     Regions mask, the Survey Depth maps and the Systematic maps.

    """
    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.
            
            intermediate_tables: a dict containing the keys:
            'good_regions' - list of dicts in which each dict must have the
             keys 'table' and 'schema'. A intersection between all this tables
             will be made based on the column pixel.
             
            'bad_regions' - list of dicts in which each dict must have the
             keys 'table' and 'schema'. The remaining pixels from the
             'good_regions' will be removed based on the bad_regions existence 
            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['good_regions', 'bad_regions']
        if not set(keys).issubset(self.intermediate_tables.keys()):
            raise ValueError("Missing intermediate table(s)")

        if len(self.intermediate_tables['good_regions']) <= 0:
            raise ValueError("At least one filter table must be added")

    def select(self):
        good_region_tables = []
        bad_region_tables = []

        for table in self.intermediate_tables['good_regions']:
            good_region_tables.append(Table(table['table'], self.dal.meta,
                                      autoload=True, schema=table['schema']))

        for table in self.intermediate_tables['bad_regions']:
            bad_region_tables.append(Table(table['table'], self.dal.meta,
                                     autoload=True, schema=table['schema']))

        stm = select([literal_column("1").label('signal'),
                      good_region_tables[0].c.pixel,
                      good_region_tables[0].c.ra,
                      good_region_tables[0].c.dec])

        # join statement
        stm_join = good_region_tables[0]
        # Inner join
        for table in good_region_tables[1:]:
            stm_join = stm_join.join(table, good_region_tables[0].c.pixel ==
                                     table.c.pixel)
        # Left Join
        for table in bad_region_tables:
            stm_join = stm_join.join(table, good_region_tables[0].c.pixel ==
                                     table.c.pixel, isouter=True)

        if len(good_region_tables) > 0 or len(bad_region_tables) > 0:
            stm = stm.select_from(stm_join)

        if len(bad_region_tables) > 0:
            for table in bad_region_tables:
                stm = stm.where(table.c.pixel == None)

        return stm


class Reduction(Operation):
    """ TODO Reduction docs
    """

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the key:
            'coadd_objects_ring' - A dict containing the keys 'table' and
            
             'schema'. This table makes the connection between maps and
             catalogs.  
            
            intermediate_tables: a dict containing the key:
            'footprint' - A dict containing the keys 'table' and
             'schema'. This table was created on the previous step

            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        if 'coadd_objects_ring' not in self.input_params:
            raise ValueError("Missing input params")

        if 'footprint' not in self.intermediate_tables:
            raise ValueError("Missing input params")

    def select(self):
        # load tables.
        t_footprint = Table(self.intermediate_tables['footprint']['table'],
                            self.dal.meta, autoload=True,
                            schema=self.intermediate_tables['footprint']['schema'])
        t_objects_ring = Table(self.input_params['coadd_objects_ring']['table'],
                               self.dal.meta, autoload=True,
                               schema=self.input_params['coadd_objects_ring']['schema'])

        # join statement
        stm_join = t_footprint
        stm_join = stm_join.join(t_objects_ring, t_footprint.c.pixel ==
                                 t_objects_ring.c.pixel)
        stm = select([t_objects_ring.c.coadd_objects_id]).\
            select_from(stm_join)
        return stm


class ZeroPoint(Operation):
    """ TODO ZeroPoint DOCS here
    """
    BANDS = ['g', 'r', 'i', 'z', 'y']
    CORRECTION_TYPES = {
        "extinction_and_slr",
        "only_extinction",
        "sfd98"
        }

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'zero_point' - A dict containing the keys 'table' and 'schema'.
            
            'coadd_objects' - A dict containing the keys 'table' and 'schema'.
            The dataset itself
            
            'correction_type' - Options: "extinction_and_slr",
            "only_extinction", "sfd98"
            
            'add_slr_shift_columns' - Boolean
            
            'add_cuts_columns' - Boolean
            
            'mag_type' - Options: 'detmodel', 'auto', 'wavgcalib'
            'aper_4', 'aper_8'
            
            'columns' - columns to apply zp correction
            """
        super().__init__(input_params, intermediate_tables)

        self.t_coadd = None
        self.t_zp = None
        self.columns = None

    def validate(self):
        keys = ['zero_point', 'coadd_objects', 'correction_type',
                'add_slr_shift_columns', 'add_cuts_columns', 'mag_type',
                'columns']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        if self.input_params['correction_type'] not in ZeroPoint.CORRECTION_TYPES:
            raise "Correction_type unvailable."

        # load tables.
        self.t_coadd = Table(self.input_params['coadd_objects']['table'], self.dal.meta,
                             autoload=True,
                             schema=self.input_params['coadd_objects']['schema']).alias('data_set')

        self.t_zp = Table(self.input_params['zero_point']['table'], self.dal.meta, autoload=True,
                          schema=self.input_params['zero_point']['schema']).alias('zero_point')

        columns_apply_zp = ZeroPoint.columns_to_apply_zp(self.input_params['columns'])
        columns_cuts = self.get_columns_from_cuts_op()
        self.columns = columns_apply_zp | columns_cuts

        corrected_columns = self.get_columns_corrected()
        slr_columns = self.get_slr_shift_corrected()

        stm_join = self.t_coadd
        stm_join = stm_join.join(self.t_zp, self.t_zp.c.coadd_objects_id ==
                                 self.t_coadd.c.coadd_objects_id)
        all_columns = [self.t_coadd.c.coadd_objects_id] + corrected_columns +\
            slr_columns
        stm = select(all_columns).select_from(stm_join)
        return stm

    def get_slr_shift_corrected(self):
        slr = []
        if self.input_params['add_slr_shift_columns']:
            for band in ZeroPoint.BANDS:
                cur_slr = ""

                if self.input_params['correction_type'] == 'extinction_and_slr':
                    col_zp_ext = getattr(self.t_zp.c, "ext_%s" % band)
                    col_zp_minus = getattr(self.t_zp.c,
                                           "zp_minus_ext_%s" % band)
                    cur_slr = - col_zp_ext + col_zp_minus
                elif self.input_params['correction_type'] == 'only_extinction':
                    col_zp_ext = getattr(self.t_zp.c, "ext_%s" % band)
                    cur_slr = - col_zp_ext
                elif self.input_params['correction_type'] == 'sfd98':
                    col_coadd_xcorr_sfd98 = getattr(self.t_coadd.c,
                                                    'xcorr_sfd98_%s' % band)
                    cur_slr = - col_coadd_xcorr_sfd98

                slr.append((cur_slr).label('slr_shift_%s' % band))
        return slr

    def get_columns_corrected(self):
        cases = []
        for column in self.columns:
            _filter = column[-1]
            col_coadd = getattr(self.t_coadd.c, column)
            cur_else = ""

            if self.input_params['correction_type'] == 'extinction_and_slr':
                col_zp_ext = getattr(self.t_zp.c, "ext_%s" % _filter)
                col_zp_minus = getattr(self.t_zp.c,
                                       "zp_minus_ext_%s" % _filter)
                cur_else = col_coadd - col_zp_ext + col_zp_minus
            elif self.input_params['correction_type'] == 'only_extinction':
                col_zp_ext = getattr(self.t_zp.c, "ext_%s" % _filter)
                cur_else = col_coadd - col_zp_ext
            elif self.input_params['correction_type'] == 'sfd98':
                col_coadd_xcorr_sfd98 = getattr(self.t_coadd.c,
                                                'xcorr_sfd98_%s' % _filter)
                cur_else = col_coadd - col_coadd_xcorr_sfd98

            cases.append(case([(col_coadd == literal_column('99'),
                         literal_column('99')), ],
                         else_=cur_else).label(column))
        return cases

    def get_columns_from_cuts_op(self):
        columns = set()
        if self.input_params['add_cuts_columns']:
            for band in Cuts.BANDS:
                columns.add(Cuts.to_magerr_column(self.input_params['mag_type'], band))
                columns.add(Cuts.to_mag_column(self.input_params['mag_type'], band))
        return columns

    @staticmethod
    def is_zero_point_column(column):
        if 'mag_' in column:
            return True
        return False

    @staticmethod
    def columns_to_apply_zp(columns):
        columns_zp = set()
        for column in columns:
            if ZeroPoint.is_zero_point_column(column):
                columns_zp.add(column)
        return columns_zp



class Cuts(Operation):
    """ TODO Cuts DOCS here
    """
    BANDS = ['g', 'r', 'i', 'z', 'y']

    MAG_TYPE = {}
    MAG_TYPE['detmodel'] = 'mag_detmodel'
    MAG_TYPE['auto'] = 'mag_auto'
    MAG_TYPE['wavgcalib'] = 'wavgcalib_mag_psf'
    MAG_TYPE['aper_4'] = 'mag_aper_4'
    MAG_TYPE['aper_8'] = 'mag_aper_8'

    MAGERR_TYPE = {}
    MAGERR_TYPE['detmodel'] = 'magerr_detmodel'
    MAGERR_TYPE['auto'] = 'magerr_auto'
    MAGERR_TYPE['wavgcalib'] = 'wavg_magerr_psf'
    MAGERR_TYPE['aper_4'] = 'magerr_aper_4'
    MAGERR_TYPE['aper_8'] = 'magerr_aper_8'

    @staticmethod
    def to_mag_column(mag_type, band):
        return ('%s_%s' % (Cuts.MAG_TYPE[mag_type], band))

    @staticmethod
    def to_magerr_column(mag_type, band):
        return ('%s_%s' % (Cuts.MAGERR_TYPE[mag_type], band))

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'coadd_objects': A dict containing the keys 'table' and 'schema'.
            The dataset itself
            
            'mag_type': Magnitude Type - options: 'detmodel', 'auto',
            'wavgcalib', 'aper_4', 'aper_8'
            
            'sextractor_flags': Sextractor reference flags. An array of any
            combination of the following FLAGS:
            0 - Clean object
            1 - The object has neighbours, bright and close enough to significantly bias the MAG_AUTO photometry, or bad pixels
            2 - The object was originally blended with another one
            4 - At least one pixel of the object is saturated (or very close to)
            8 - The object is truncated (too close to an image boundary)
            16 - Object aperture data are incomplete or corrupted
            32 - Object isophotal data are incomplete or corrupted
            64 - A memory overflow occurred during deblending
            128 - A memory overflow occurred during extraction
            
            'sextractor_bands': Sextractor reference bands - Any combination
             of the bands: 'g','r','i','z','Y'.
             
            'additional_cuts': Additional cuts - Any combination of the
            following parameters:
            remove_bbj - Remove Bright Blue Junks 
            spreaderr_model - Remove objects where the spreadmodel fit failed 
            bad_astronomic_colors - Remove objects with bad astrometric colors
             
            'niter_model': Niter Model - Any combination of the bands:
            'g','r','i','z','Y'.
            
            'sn_cuts': Signal-to-noise cuts - An dict associating band and
            value. e.g.: TODO
            
            'magnitude_limit': Magnitude limit cuts - An dict associating band
            and value. e.g.: {'i': '22'}
            
            'bright_magnitude': Bright magnitude cuts - An dict associating
            band and value. e.g.: {'i': '15'}
            
            'color_cuts': Color cuts - An dict associating -'gr', 'ri', 'iz', 
            'zy'- and -min, max- values. e.g.: {'gr':[-2, 4], 'ri':[-2, 4]}

            intermediate_tables: a dict containing the key:
            'zero_point' - A dict containing the keys 'table' and
             'schema'. This table was created on the previous step
            'reduction' - A dict containing the keys 'table' and
             'schema'. This table was created on the previous step

            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['coadd_objects', 'mag_type', 'sextractor_flags',
                'sextractor_bands', 'additional_cuts', 'niter_model',
                'sn_cuts', 'magnitude_limit', 'bright_magnitude',
                'color_cuts']

        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        t_reduction = Table(self.intermediate_tables['reduction']['table'],
                            self.dal.meta, autoload=True,
                            schema=self.intermediate_tables['reduction']['schema']).alias('reduction')
        t_coadd = Table(self.input_params['coadd_objects']['table'], self.dal.meta,
                        autoload=True,
                        schema=self.input_params['coadd_objects']['schema']).alias('coadd_objects')

        # join statement
        stm_join = t_reduction
        stm_join = stm_join.join(t_coadd, t_reduction.c.coadd_objects_id ==
                                 t_coadd.c.coadd_objects_id)

        _where = []

        # cuts involving only coadd_objects_columns
        # sextractor flags
        if 'sextractor_bands' in self.input_params and\
           'sextractor_flags' in self.input_params:
            # combine_flags
            queries = []
            sum_flags = sum(self.input_params['sextractor_flags'])
            for band in self.input_params['sextractor_bands']:
                query = []
                col = getattr(t_coadd.c, 'flags_%s' % band)
                if 0 in self.input_params['sextractor_flags']:
                    query.append(col == literal_column('0'))
                if sum_flags > 0:
                    and_op = se.BitwiseAnd(
                                   col,
                                   literal_column(str(sum_flags)))
                    query.append((and_op) > literal_column('0'))
                queries.append(or_(*query))
            _where.append(and_(*queries))

        # bbj
        if 'remove_bbj' in self.input_params['additional_cuts']:
            _where.append(or_(
                            t_coadd.c.nepochs_g > literal_column('0'),
                            t_coadd.c.magerr_auto_g > literal_column('0.05'),
                            t_coadd.c.mag_model_i - t_coadd.c.mag_auto_i >
                            literal_column('-0.4')
                            ))

        # niter model
        if 'niter_model' in self.input_params:
            tmp = []
            for band in self.input_params['niter_model']:
                col = getattr(t_coadd.c, 'niter_model_%s' % band)
                tmp.append(col > literal_column('0'))
            _where.append(and_(*tmp))

        # spreaderr model
        if 'spreaderr_model' in self.input_params['additional_cuts']:
            tmp = []
            for band in Cuts.BANDS:
                col = getattr(t_coadd.c, 'spreaderr_model_%s' % band)
                tmp.append(col > literal_column('0'))
            _where.append(and_(*tmp))

        # bad astronomic color
        if 'bad_astronomic_colors' in self.input_params['additional_cuts']:
            _where.append(or_(
                    and_(
                        func.abs(t_coadd.c.alphawin_j2000_g -
                                 t_coadd.c.alphawin_j2000_i) <
                        literal_column('0.0003'),
                        func.abs(t_coadd.c.deltawin_j2000_g -
                                 t_coadd.c.deltawin_j2000_i) <
                        literal_column('0.0003')
                    ),
                    t_coadd.c.magerr_auto_g > literal_column('0.05')

                ))

        # ra, dec
        if 'dec_min' in self.input_params:
            _where.append(t_coadd.c.dec >= self.input_params['dec_min'])
        if 'dec_max' in self.input_params:
            _where.append(t_coadd.c.dec <= self.input_params['dec_max'])
        if 'ra_min' in self.input_params:
            _where.append(t_coadd.c.ra >= self.input_params['ra_min'])
        if 'ra_max' in self.input_params:
            _where.append(t_coadd.c.ra <= self.input_params['ra_max'])

        t_cur = t_coadd
        if 'zero_point' in self.intermediate_tables:
            t_cur = Table(self.intermediate_tables['zero_point']['table'], self.dal.meta,
                          autoload=True, schema=self.intermediate_tables['zero_point']['schema']).\
                          alias('zero_point')
            stm_join = stm_join.join(t_cur, t_reduction.c.coadd_objects_id ==
                                     t_cur.c.coadd_objects_id)

        # signal to noise cuts
        if 'sn_cuts' in self.input_params:
            tmp = []
            for element in self.input_params['sn_cuts'].items():
                band, value = element
                db_col = Cuts.to_magerr_column(self.input_params['mag_type'], band)
                col = getattr(t_cur.c, db_col)
                tmp.append(and_(
                        col > literal_column('0'),
                        literal_column('1.086')/col >
                        literal_column(str(value))
                    ))
            _where.append(and_(*tmp))

        # magnitude limit
        if 'magnitude_limit' in self.input_params:
            tmp = []
            for element in self.input_params['magnitude_limit'].items():
                band, value = element
                db_col = Cuts.to_mag_column(self.input_params['mag_type'], band)
                col = getattr(t_cur.c, db_col)
                tmp.append(col < literal_column(str(value)))
            _where.append(and_(*tmp))

        # bright magnitude limit
        if 'bright_magnitude' in self.input_params:
            tmp = []
            for element in self.input_params['bright_magnitude'].items():
                band, value = element
                db_col = Cuts.to_mag_column(self.input_params['mag_type'], band)
                col = getattr(t_cur.c, db_col)
                tmp.append(col > literal_column(str(value)))
            _where.append(and_(*tmp))

        # color cuts
        if 'color_cuts' in self.input_params:
            tmp = []
            for element in self.input_params['color_cuts'].items():
                band, value = element
                db_col_max = Cuts.to_mag_column(self.input_params['mag_type'], band[0])
                db_col_min = Cuts.to_mag_column(self.input_params['mag_type'], band[1])
                col_max = getattr(t_cur.c, db_col_max)
                col_min = getattr(t_cur.c, db_col_min)
                tmp.append(between(literal_column(str(col_max - col_min)),
                                   literal_column(str(value[0])),
                                   literal_column(str(value[1]))))
            _where.append(and_(*tmp))

        stm = select([t_coadd.c.coadd_objects_id]).\
            select_from(stm_join).where(and_(*_where))

        return stm


class Bitmask(Operation):
    """ TODO Bitmask DOCS here
    """
    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'molygon_coadds' - A dict containing the keys 'table' and 'schema'.
            
            'molygon' - A dict containing the keys 'table' and 'schema'.
            
            'mangle_bitmask' - Options:  Any combination of the bands:
            'g','r','i','z','Y'.
            
            intermediate_tables: a dict containing the key:
            'last_table' - A dict containing the keys 'table' and
             'schema'. The table created on the previous step
            
            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['molygon_coadds', 'molygon', 'mangle_bitmask']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        sub_op = list(self.intermediate_tables.values())[0]

        # load tables.
        t_sub_op = Table(self.intermediate_tables['last_table']['table'], self.dal.meta, autoload=True,
                         schema=self.intermediate_tables['last_table']['schema'])
        _where = []

        # bitmask
        alias_table = None
        t_coadd_molygon = Table(self.input_params['molygon_coadds']['table'],
                                self.dal.meta, autoload=True,
                                schema=self.input_params['molygon_coadds']['schema']).alias('molygon_coadds')
        t_molygon = Table(self.input_params['molygon']['table'], self.dal.meta,
                          autoload=True, schema=self.input_params['molygon']['schema']).alias('molygon')

        stm_join = t_sub_op
        stm_join = stm_join.join(t_coadd_molygon,
                                 t_sub_op.c.coadd_objects_id ==
                                 t_coadd_molygon.c.coadd_objects_id)

        for band in self.input_params['mangle_bitmask']:
            # give the str column and retrieve the attribute.
            alias_table = t_molygon.alias('molygon_%s' % band)
            col = getattr(t_coadd_molygon.c, 'molygon_id_%s' % band)
            stm_join = stm_join.join(alias_table,
                                     col == alias_table.c.id)
        _where.append(alias_table.c.hole_bitmask != literal_column('1'))

        stm = select([t_sub_op.c.coadd_objects_id]).\
            select_from(stm_join).where(and_(*_where))

        return stm


class ObjectSelection(Operation):
    """ TODO ObjectSelection DOCS here
    """

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'coadd_objects' - A dict containing the keys 'table' and 'schema'.

            'columns_data_set' - List of columns present on the dataset table.

            'columns_zero_point' - Options:  Any combination of the bands:
            'g','r','i','z','Y'.

            intermediate_tables: a dict containing the key:
            'last_table' - A dict containing the keys 'table' and
             'schema'. The table created on the previous step

            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['coadd_objects', 'columns_data_set']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        sub_op = list(self.intermediate_tables.values())[0]

        # load tables.
        t_table = Table(self.intermediate_tables['last_table']['table'], self.dal.meta, autoload=True,
                        schema=self.intermediate_tables['last_table']['schema'])
        columns = [t_table.c.coadd_objects_id]

        stm_join = t_table
        if self.input_params['columns_data_set']:
            t_coadd = Table(self.input_params['coadd_objects']['table'], self.dal.meta,
                            autoload=True,
                            schema=self.input_params['coadd_objects']['schema']).alias('coadd_objects')
            stm_join = stm_join.join(t_coadd,
                                     t_table.c.coadd_objects_id ==
                                     t_coadd.c.coadd_objects_id)
            for column in self.input_params['columns_data_set']:
                columns.append(getattr(t_coadd.c, column))

        if self.input_params['columns_zero_point']:
            t_zp = Table(self.intermediate_tables['zero_point']['table'], self.dal.meta, autoload=True,
                         schema=self.intermediate_tables['zero_point']['schema']).alias('zero_point')
            stm_join = stm_join.join(t_zp,
                                     t_table.c.coadd_objects_id ==
                                     t_zp.c.coadd_objects_id)
            for column in self.input_params['columns_zero_point']:
                columns.append(getattr(t_zp.c, column))

        stm = select(columns).select_from(stm_join)
        return stm


class SgSeparation(Operation):
    """ TODO SgSeparation DOCS here
    """

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'select_starts' -

            'select_galaxies' -

            'reference_band' - A dict containing the keys 'table' and
             'schema'.
             
            'sg_tables' - A list of dicts containing the keys 'table' and
             'schema'.

            intermediate_tables: a dict containing the key:
            'last_table' - A dict containing the keys 'table' and
             'schema'. The table created on the previous step

            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['select_starts', 'select_galaxies', 'reference_band',
                'sg_tables']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        # load tables.
        t_obj_selection = Table(self.intermediate_tables['last_table']['table'],
                                self.dal.meta, autoload=True,
                                schema=self.intermediate_tables['last_table']['schema'])
        t_sg = []
        for table in self.input_params['sg_tables']:
            t_sg.append(Table(table['table'], self.dal.meta, autoload=True,
                              schema=table['schema']))

        _where = []
        # join statement
        stm_join = t_obj_selection
        for table in t_sg:
            stm_join = stm_join.join(
                table, t_obj_selection.c.coadd_objects_id ==
                table.c.coadd_objects_id)
            col = getattr(table.c, '%s' % self.input_params['reference_band'])
            if self.input_params['select_galaxies']:
                _where.append(col == literal_column('0'))
            if self.input_params['select_starts']:
                _where.append(col == literal_column('1'))

        stm = select([t_obj_selection.c.coadd_objects_id]).\
            select_from(stm_join).where(and_(*_where))

        return stm

    #REVIEW - case when both star and galaxies are selected.
    def select_all_objects(self, input_params):
        if input_params['select_starts'] and input_params['select_galaxies']:
            return True
        return False


class PhotoZ(Operation):
    """ TODO PhotoZ DOCS here
    """

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'zmin' -
             
            'zmax' -
            
            'pz_tables' - A list of dicts containing the keys 'table' and
             'schema'.

            intermediate_tables: a dict containing the key:
            'last_table' - A dict containing the keys 'table' and
             'schema'. The table created on the previous step

            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['zmin', 'zmax', 'pz_tables']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        # load tables.
        t_sub_op = Table(self.intermediate_tables['last_table']['table'],
                         self.dal.meta, autoload=True,
                         schema=self.intermediate_tables['last_table']['schema'])
        t_pz = []
        for table in self.input_params['pz_tables']:
            t_pz.append(Table(table['table'], self.dal.meta, autoload=True,
                              schema=table['schema']))

        _where = []
        # join statement
        stm_join = t_sub_op
        for table in t_pz:
            stm_join = stm_join.join(
                table, t_sub_op.c.coadd_objects_id ==
                table.c.coadd_objects_id)
            _where.append(and_(table.c.z_best >
                               literal_column(str(self.input_params['zmin'])),
                               table.c.z_best <
                               literal_column(str(self.input_params['zmax']))))

        stm = select([t_sub_op.c.coadd_objects_id]).\
            select_from(stm_join).where(and_(*_where))

        return stm


class GalaxyProperties(Operation):
    """ TODO GalaxyProperties DOCS here
    """

    def __init__(self, input_params, intermediate_tables):
        """
        Args:
            dal: database connection class.

            input_params: a dict containing the keys:
            'columns' - Galaxy properties columns
            
            'gp_tables' - A list of dicts containing the keys 'table' and
             'schema'.

            intermediate_tables: a dict containing the key:
            'last_table' - A dict containing the keys 'table' and
             'schema'. The table created on the previous step

            """
        super().__init__(input_params, intermediate_tables)

    def validate(self):
        keys = ['columns', 'gp_tables']
        if not set(keys).issubset(self.input_params.keys()):
            raise ValueError("Missing input_param(s)")

    def select(self):
        # load tables.
        t_sub_op = Table(self.intermediate_tables['last_table']['table'], self.dal.meta, autoload=True,
                         schema=self.intermediate_tables['last_table']['schema'])
        t_gp = []
        for table in self.input_params['gp_tables']:
            t_gp.append(Table(table['table'], self.dal.meta, autoload=True,
                              schema=table['schema']))

        # join statement
        stm_join = t_sub_op
        for table in t_gp:
            stm_join = stm_join.join(
                table, t_sub_op.c.coadd_objects_id == table.c.coadd_objects_id)

        # REVIEW - only the first column is being considered for
        columns = [t_sub_op.c.coadd_objects_id]
        for column in self.input_params['columns']:
            columns.append(getattr(t_gp[0].c, column))

        stm = select(columns).select_from(stm_join)
        return stm
