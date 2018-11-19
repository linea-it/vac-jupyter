fracdet_g = (
    'CREATE TEMP TABLE fracdet_g AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031864 a)')

fracdet_r = (
    'CREATE TEMP TABLE fracdet_r AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031865 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 0.9)')

fracdet_i = (
    'CREATE TEMP TABLE fracdet_i AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031851 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 0.9)')

fracdet_z = (
    'CREATE TEMP TABLE fracdet_z AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031866 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 0.9)')

fracdet_Y = (
    'CREATE TEMP TABLE fracdet_Y AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031867 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 0.9)')

fracdet_all = (
    'CREATE TEMP TABLE fracdet_all AS '
    '(SELECT a.pixel, 1 as signal, a.ra, a.dec FROM fracdet_g a '
    'INNER JOIN fracdet_r b ON a.pixel= b.pixel '
    'INNER JOIN fracdet_i c ON a.pixel= c.pixel '
    'INNER JOIN fracdet_z d ON a.pixel= d.pixel '
    'INNER JOIN fracdet_y e ON a.pixel= e.pixel)')

depth_map_auto_i = (
    'CREATE TEMP TABLE depth_map_auto_i AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec FROM upload.fits_10031873 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 21.0)')

depth_map_auto_r = (
    'CREATE TEMP TABLE depth_map_auto_r AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec FROM upload.fits_10031872 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'from y3_gold.y3_gold_2_2_small) and signal >= 23.0)')

depth_map_auto_all = (
    'CREATE TEMP TABLE depth_map_auto_all AS '
    '(SELECT a.pixel, 1 as signal, a.ra, a.dec FROM depth_map_auto_i a '
    'INNER JOIN depth_map_auto_r b ON a.pixel= b.pixel)'
)

n_images_g = (
    'CREATE TEMP TABLE n_images_g AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031878 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 1)')

n_images_r = (
    'CREATE TEMP TABLE n_images_r AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031879 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 1)')

n_images_i = (
    'CREATE TEMP TABLE n_images_i AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec '
    'FROM upload.fits_10031897 a '
    'WHERE pixel in (SELECT distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and signal >= 1)')

n_images_all = (
    'CREATE TEMP TABLE n_images_all AS '
    '(SELECT a.pixel, 1 as signal, a.ra, a.dec '
    'FROM n_images_g a '
    'INNER JOIN n_images_r b ON a.pixel= b.pixel '
    'INNER JOIN n_images_i c ON a.pixel= c.pixel)'
)

bad_regions = (
    'CREATE TEMP TABLE bad_regions AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec FROM upload.fits_10031895 a '
    'WHERE pixel in (select distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and (CAST(signal AS INTEGER) & 3) > 0)'
)

foreground = (
    'CREATE TEMP TABLE foreground AS '
    '(SELECT a.pixel, a.signal, a.ra, a.dec FROM upload.fits_10031896 a '
    'WHERE pixel in (select distinct(hpix_4096) '
    'FROM y3_gold.y3_gold_2_2_small) and (CAST(signal AS INTEGER) & 127) > 0)')

footprint = (
    'CREATE TEMP TABLE footprint AS '
    '(SELECT a.pixel, 1 AS signal, a.ra, a.dec, '
    'upload.fits_10031864.signal AS detfrac_g, '
    'upload.fits_10031865.signal AS detfrac_r, '
    'upload.fits_10031851.signal AS detfrac_i, '
    'upload.fits_10031866.signal AS detfrac_z, '
    'upload.fits_10031867.signal AS detfrac_Y '
    'FROM fracdet_all a '
    'INNER JOIN depth_map_auto_i b ON a.pixel= b.pixel '
    'INNER JOIN n_images_all c ON a.pixel= c.pixel '  
    'LEFT JOIN bad_regions d ON a.pixel = d.pixel '
    'LEFT JOIN foreground e ON a.pixel = e.pixel '
    'LEFT JOIN upload.fits_10031864 ON a.pixel = upload.fits_10031864.pixel '
    'LEFT JOIN upload.fits_10031865 ON a.pixel = upload.fits_10031865.pixel '
    'LEFT JOIN upload.fits_10031851 ON a.pixel = upload.fits_10031851.pixel '
    'LEFT JOIN upload.fits_10031866 ON a.pixel = upload.fits_10031866.pixel '
    'LEFT JOIN upload.fits_10031867 ON a.pixel = upload.fits_10031867.pixel ' 
    'WHERE d.pixel IS NULL and e.pixel IS NULL)'
)

reduction = (
    'CREATE TEMP TABLE reduction AS (SELECT a.mag_auto_g AS mag_g, '
    'a.mag_auto_r AS mag_r, '
    'a.mag_auto_i AS mag_i, '
    'a.mag_auto_z AS mag_z, '
    'a.mag_auto_y AS mag_y, '
    'a.magerr_auto_g AS magerr_g, '
    'a.magerr_auto_r AS magerr_r, '
    'a.magerr_auto_i AS magerr_i, '
    'a.magerr_auto_y AS magerr_y, '
    'a.magerr_auto_z AS magerr_z, '
    'a.coadd_objects_id, '
    'a.ra, '
    'a.dec, '
    'a.hpix_4096 '
    'FROM y3_gold.y3_gold_2_2_des_small a '
    'INNER JOIN footprint b ON a.hpix_4096 = b.pixel '
    'WHERE (NOT(((CAST(flags_i AS INTEGER) & 255) > 0))) '
    'AND (extended_class_mash_mof = 2 '
    'OR extended_class_mash_mof = 3))'
)

object_selection = (
    'CREATE TEMP TABLE object_selection AS (SELECT * from reduction '
    'WHERE (mag_i > 15) '
    'AND (mag_i < 23.0) '
    'AND (mag_g - mag_r BETWEEN -2.0 AND 4.0 '
    'AND mag_r - mag_i BETWEEN -2.0 AND 4.0 '
    'AND mag_i - mag_z BETWEEN -2.0 AND 4.0 '
    'AND mag_z - mag_Y BETWEEN -2.0 AND 4.0))'
)

footprint_checked = (
    'CREATE TEMP TABLE footprint_checked AS ( SELECT a.* '
    'FROM footprint a '
    'WHERE pixel IN '
    '(SELECT distinct(hpix_4096) '
    'FROM object_selection))'
)

def create_btree_index(table_name, column='pixel'):
    return 'CREATE INDEX %(tn)s_%(col)s ON %(tn)s ' \
           'USING btree(%(col)s)' % ({'tn': table_name, 'col': column})