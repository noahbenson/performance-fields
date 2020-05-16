####################################################################################################
# Visual Performance Fields Dataset Description
# This file implements a neuropythy dataset object fo the visual performance
# fields project by Benson, Kupers, Carrasco, and Winawer.
# By Noah C. Benson

import neuropythy, pimms

neuropythy.config.declare('visual_performance_fields_path',
                          environ_name='VISUAL_PERFORMANCE_FIELDS_PATH',
                          default_value=None, filter=neuropythy.datasets.hcp.to_nonempty_path)
@pimms.immutable
class VisualPerformanceFieldsDataset(neuropythy.datasets.HCPMetaDataset):
    '''
    The VisualPerformanceFieldsDataset class declares the basic data structures and the functions
    for generating, loading, or preprocessing the visual performance fields dataset associated
    with the paper by Benson, Kupers, Carrasco, and Winawer (2020).
    
    This dataset contains the following members:
      * subject_list is a list of the subject IDs for each HCP subject analyzed in this notebook.
      * metadata_table is a pandas dataframe of all the meta-data (i.e., the HCP behavioral data)
        for each subject.
      * gender is a dictionary whose keys are subject IDs and whose values are either 'M' or 'F'.
      * agegroup is a dictionary of age-groups for each subject; the age-group value is the age in
        the middleof the age range for a particular subject. E.g., if  subject u is in the 26-30
        age-group, then agegroup[u] will be 28.
      * inferred_maps is a dictionary mapping subject IDs to retinotopic mapping properties as
        inferred by Bayesian inference (Benson and Winawer, 2018, DOI:10.7554/eLife.40224).
      * cleaned_maps is a equivalent to inferred_maps but contains inferred retinotopic
        parameters after a small amount of smoothing designed to even-out improbably cortical
        magnification values. The value cleaned_maps[sid][h][p][k] is the property p ('polar_angle',
        'eccentricity', 'radius', or 'visual_area') for vertex k*of the hemisphere h ('lh' or 'rh')
        of subject sid.
      * v1_distance is a dictionary mapping subject IDs to properties that specify the distance from
        each (ventral/dorsal) V1 boundary. The value v1_distance[sid][h][vd][k] is the
        midgray-surface distance from vertex k to the vd ('ventral' or 'dorsal') V1 boundary of
        hemisphere h ('lh' or 'rh') of subject sid.
      * inferred_table is a dataframe summarizing the inferred retinotopic maps and the various data
        related to them (cleaned data and distances).
      * subjects is a dictionary of subject objects, each of which contain, as properties or
        meta-data, the various data described above.
    '''
    
    # Constants / Magic Numbers ####################################################################
    roi_angles = (10, 20, 30, 40, 50)
    roi_angles_fine = (5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    # The eccentricity ranges we look at:
    roi_eccens = ((0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7))

    def __init__(self, url=None, cache_directory=Ellipsis,
                 meta_data=None, create_directories=True, create_mode=0o755):
        cdir = cache_directory
        if cdir is Ellipsis:
            cdir = neuropythy.config['visual_performance_fields_path']
        neuropythy.datasets.HCPMetaDataset.__init__(self, name='visual_performance_fields',
                                                    cache_directory=cdir,
                                                    meta_data=meta_data, create_mode=create_mode,
                                                    create_directories=create_directories)
        self.url = url

    # Administrative Data ##########################################################################
    # Things that are required to load or organize the data but that aren't the data themselves.
    @pimms.param
    def url(url):
        '''
        url is the URL from which the performance-fields data is loaded.
        '''
        if url is None or url is Ellipsis: return 'osf://5gprz/'
        if not pimms.is_str(u): raise ValueError('url must be a string')
        return u
    @pimms.value
    def pseudo_path(url, cache_directory):
        '''
        pseudo_path is the neuropythy psueod-path object responsible for loading the OSF data
        related to the visual performance fields dataset.
        '''
        from neuropythy.util import pseudo_path
        return pseudo_path(url, cache_path=cache_directory).persist()
    inferred_map_files = {'polar_angle':  '%s.%s.inf-MSMAll_angle.native32k.mgz',
                          'eccentricity': '%s.%s.inf-MSMAll_eccen.native32k.mgz',
                          'radius':       '%s.%s.inf-MSMAll_sigma.native32k.mgz',
                          'visual_area':  '%s.%s.inf-MSMAll_varea.native32k.mgz'}
    @pimms.value
    def hcp_data():
        '''
        hcp_data is the HCP Retinotopy dataset that is used in conjunction with the Visual
        Performance Fields dataset.
        '''
        return neuropythy.data['hcp_retinotopy']
    @pimms.value
    def subject_list(hcp_data):
        '''
        subjcet_list is the list of subjects used in the dataset.
        '''
        sids = [sid for sid in hcp_data.subject_ids if sid < 999990]
        return tuple(sids)

    # Supporting Data ##############################################################################
    # Data that is critical for calculating the wedge ROIs but that aren't themselves the data used
    # to make the paper plots.
    @pimms.value
    def inferred_maps(pseudo_path, subject_list):
        '''
        inferred_maps is a nested-dictionary structure containing the retinotopic maps inferred by
        using Bayesian inference on the retinotopic maps of the subjects in the HCP 7T Retinotopy
        Dataset.
        '''
        import os, six
        from neuropythy.util import curry
        from neuropythy import load
        inffiles = VisualPerformanceFieldsDataset.inferred_map_files
        def _load_infmaps(sid,h,patt):
            flnm = pseudo_path.local_path('inferred_maps', patt % (sid,h))
            return load(flnm)
        return pimms.persist(
            {sid: {h: pimms.lmap({('inf_'+k): curry(_load_infmaps, sid, h, v)
                                  for (k,v) in six.iteritems(inffiles)})
                   for h in ['lh','rh']}
             for sid in subject_list})
    # How to calculate the cleaned maps
    @staticmethod
    def _generate_clean_maps(sid, h, infmaps):
        import neuropythy as ny
        import numpy as np
        sub = ny.hcp_subject(sid)
        hemi = sub.hemis[h]
        # we use the inf-lowres-prf_* retinotopy data
        irdat = infmaps[sid][h]
        # get x/y and labels
        lbl = irdat['inf_visual_area']
        (iang,iecc) = ny.as_retinotopy(irdat, 'visual')
        # get the mask
        mask = np.where((iecc <= 7) & ((lbl == 1) | (lbl == 2)))[0]
        # calculate the clean retinotopy
        (ang, ecc) = ny.vision.clean_retinotopy(
            hemi, 'prf_', visual_area=None, steps=[500]*5,
            fieldsign_knob=6, rt_knob=None,
            mask=mask, rounds=5, average=None, jitter=None)
        res = {'clean_polar_angle':ang, 'clean_eccentricity':ecc}
        return res
    @pimms.value
    def cleaned_maps(pseudo_path, subject_list):
        '''
        cleaned_maps is a nested-dictionary structure containing cleaned versions of the
        measured retinotopic maps of the subjects in the HCP 7T Retinotopy Dataset.
        '''
        import os, six
        from neuropythy.util import curry
        from neuropythy import load
        def _load_cleanmaps(sid,h):
            flnm = pseudo_path.local_path('clean_maps', '%s_%s.mgz' % (sid,h))
            (ang,ecc) = load(flnm)
            return pimms.persist({'clean_polar_angle': ang, 'clean_eccentricity': ecc})
        return pimms.persist({sid: pimms.lmap({h: curry(_load_cleanmaps, sid, h)
                                               for h in ['lh','rh']})
                              for sid in subject_list})
    @staticmethod
    def _generate_boundary_distances(sid, h, infmaps):
        import neuropythy as ny
        import numpy as np
        import six
        sub = neuropythy.hcp_subject(sid)
        hemi = sub.hemis[h]
        # we use the inferred retinotopy data loaded above
        rdat = ny.retinotopy_data(infmaps[sid][h], 'inf_')
        # get x/y and labels
        (x0,y0) = ny.as_retinotopy(rdat, 'geographical')
        lbl = rdat['visual_area']
        (ang,ecc) = ny.as_retinotopy(rdat, 'visual')
        # get V1/V2 lines:
        mask = (ecc <= 7) & (lbl > 0)
        v1 = np.where(mask & (lbl == 1))[0]
        v2 = np.where(mask & (lbl == 2))[0]
        v12 = np.union1d(v1, v2)
        # invert x for v2
        x = np.array(x0)
        y = y0
        x[v2] = -x[v2]
        # make a visual field mesh 
        vmesh = hemi.surfaces['midgray'].submesh(v12)
        vmesh = vmesh.copy(coordinates=(x[vmesh.labels],y[vmesh.labels]))
        # draw a line from origin to the end of the stimulus; see what intersects it
        addrs = {}
        paths = {}
        for (ept, nm) in zip([(0,7), (0,-7), (7,0)], ['ventral', 'dorsal', 'horizontal']):
            seg = np.array([(0,0), ept])
            ipts = ny.geometry.segment_intersection_2D(seg, vmesh.edge_coordinates)
            ipts = np.transpose(ipts)
            ipts = ipts[np.isfinite(ipts[:,0])]
            # these are the points to address
            dists = np.dot(ipts, seg[1] - seg[0])
            ipts = ipts[np.argsort(dists)]
            # possibly the first/last few are not in the mesh; just ignore these
            addr = vmesh.address(ipts)
            ii = np.isfinite(addr['coordinates'][0])
            addr = {k:v[:,ii] for (k,v) in six.iteritems(addr)}
            addrs[nm] = addr
            paths[nm] = ny.geometry.Path(hemi, addr)
        # use these paths to calculate distances
        eds = {k:p.estimated_distances['midgray'] for (k,p) in six.iteritems(paths)}
        return eds
    def generate_boundary_distances(self, sid, h):
        '''
        generate_boundary_distances(sid, h) is a method of the visual performance field dataset that
          recalculates the boundary distances for subject sid, hemisphere h, using the data found
          in the subject's retinotopic maps and inferred maps.
        '''
        f = VisualPerformanceFieldsDataset._generate_boundary_distances
        return f(sif, h, self.inferred_maps)
    @pimms.value
    def boundary_distances(pseudo_path, subject_list, inferred_maps):
        '''
        boundary_distances is a nested-dictionary structure containing distances between
        each vertex and a V1 boundary. If x is boundar_distances[sid][h][b][k] then x is
        the distance between the k'th vertex and boundary b ("ventral", "dorsal", or
        "horizontal") in the h hemisphere ("lh" or "rh") of the subject with ID sid.
        '''
        import os, six
        from neuropythy.util import curry
        from neuropythy import load
        def _load_distances(sid, h):
            flnm = pseudo_path.local_path('distances', '%s_%s.mgz' % (sid,h))
            (v,d,h) = load(flnm).T
            return pimms.persist({'ventral': v, 'dorsal': d, 'horizontal': h})
        return pimms.persist({sid: pimms.lmap({h: curry(_load_distances, sid, h)
                                               for h in ['lh','rh']})
                              for sid in subject_list})
    @staticmethod
    def _generate_summary_table(infmaps, cleanmaps, bdists, subject_list):
        import neuropythy as ny
        import six, numpy as np
        cols = ['sid','hemi','prf_polar_angle', 'prf_eccentricity', 'prf_variance_explained',
                'prf_radius', 'visual_area', 'inf_polar_angle', 'inf_eccentricity',
                'ventral', 'dorsal','ventral_distance', 'dorsal_distance']
        df = {k:[] for k in cols}
        for sid in data.subject_list:
            sub = ny.hcp_subject(sid)
            for h in ('lh','rh'):
                hemi  = sub.hemis[h]
                dists = bdists[sid][h]
                if dists is None or any(v is None for v in six.itervalues(dists)): continue
                rmaps = ny.retinotopy_data(hemi, 'prf_')
                infmap = infmaps[sid][h]
                lbls  = np.array(infmap['inf_visual_area'])
                iang  = np.array(infmap['inf_polar_angle'])
                iecc  = np.array(infmap['inf_eccentricity'])
                mask  = np.where(np.isin(lbls, [1,2]) & (iecc <= 7))[0]
                vnt   = (np.abs(iang) < 90)
                drs   = ~vnt
                n = len(mask)
                dat = {'sid': np.full(n, sid), 'hemi':np.full(n, h),
                       'prf_polar_angle': rmaps['polar_angle'][mask],
                       'prf_eccentricity': rmaps['eccentricity'][mask],
                       'prf_radius': rmaps['radius'][mask],
                       'prf_variance_explained': rmaps['variance_explained'][mask],
                       'visual_area':lbls[mask],
                       'inf_polar_angle':np.abs(iang[mask]),
                       'inf_eccentricity':iecc[mask],
                       'ventral':vnt[mask],
                       'dorsal':drs[mask],
                       'ventral_distance':dists['V1_ventral'][mask],
                       'dorsal_distance':dists['V1_dorsal'][mask]}
                for (k,v) in six.iteritems(dat): 
                    df[k].append(v)
        df = {k:np.concatenate(v) for (k,v) in six.iteritems(df)}
        return ny.to_dataframe(df)
    def generate_summary_table(self):
        '''
        generate_summary_table() is a method for the visual performance fields dataset that
        regenerates the summary table from the inferred maps and boundary distances in the
        dataset.
        '''
        f = VisualPerformanceFieldsDataset._generate_summary_table
        return f(self.inferred_maps, self.cleaned_maps, self.boundary_distances, self.subject_list)
    @pimms.value
    def summary_table(subject_list, pseudo_path):
        '''
        summary_table is a dataframe that summarizes all of the data employed in the project for the
        calculation of the various wedge ROIs. The table includes data for every vertex in the V1/V2
        ROIs (between 0 and 7 degrees of eccentricity) of both hemispheres of all subjects.
        '''
        from neuropythy import load
        flnm = pseudo_path.local_path('summary_table.csv')
        return load(flnm)
    @staticmethod
    def _generate_hemi(sid, h, infmaps, cleanmaps, bdists):
        import neuropythy as ny, six
        hem = ny.hcp_subject(sid).hemis[h]
        cl = cleanmaps[sid][h]
        bi = infmaps[sid][h]
        ds = bdists[sid][h]
        ps = dict(cl)
        for (k,v) in six.iteritems(ds): ps[k + '_distance'] = v
        for (k,v) in six.iteritems(bi): ps[k] = v
        return hem.with_prop(ps)
    @staticmethod
    def _generate_subject(sid, infmaps, cleanmaps, bdists):
        import neuropythy as ny, six
        lh = VisualPerformanceFieldsDataset._generate_hemi(sid, 'lh', infmaps, cleanmaps, bdists)
        rh = VisualPerformanceFieldsDataset._generate_hemi(sid, 'rh', infmaps, cleanmaps, bdists)
        return ny.hcp_subject(sid).with_hemi(lh=lh, rh=rh)
    @pimms.value
    def subjects(inferred_maps, cleaned_maps, boundary_distances, subject_list):
        '''
        subjects is a dictionary of subject objects for all subjects used in the visual performance
        fields dataset. All subject objects in the subejcts dict include property data on the native
        hemispheres for inferred and cleaned retinotopic maps and for V1 boundary distances.
        '''
        from neuropythy.util import curry
        f = VisualPerformanceFieldsDataset._generate_subject
        return pimms.lmap({sid: curry(f, sid, inferred_maps, cleaned_maps, boundary_distances)
                           for sid in subject_list})

    # ROI-calculation Functions ####################################################################
    # These are methods that calculate the distance-based ROIs.
    @staticmethod
    def _generate_DROIs(sub, h, dist_prop, val_prop, ref_val, val_diff,
                        masks=None, distance_masks=None, inv_prop=None):
        import neuropythy as ny, numpy as np, six
        hem = sub.hemis[h]
        # Get the distance:
        d0 = np.abs(hem.property(dist_prop)) # (abs makes a copy)
        if inv_prop is not None:
            x = np.abs(hem.property(inv_prop))
            masks.append(np.where(d0 < x)[0])
        # Get the property:
        v0 = hem.property(val_prop)
        if val_prop == 'prf_polar_angle' and h == 'rh': v0 = -v0
        v  = np.abs(v0 - ref_val) - val_diff
        # Apply the distance masks; the masks are for deciding which vertices get
        # returned as part of the ROI while the distance masks are for calculating
        # the distance boundary; for the distance boundary we actually need both
        ii = hem.indices
        for mask in ([] if masks is None else masks):
            ii = np.intersect1d(ii, hem.mask(mask))
        mask_ii = ii
        for mask in ([] if distance_masks is None else distance_masks):
            ii = np.intersect1d(ii, hem.mask(mask))
        # clear out values that aren't in the distance mask; we clear them out of
        # v (the value on which we're searching over for the distance)
        ii = np.setdiff1d(hem.indices, ii)
        v[ii] = np.nan
        # Find edges that straddle the val_diff from the ref_val
        (a,b) = hem.tess.indexed_edges
        v = np.sign(v[a] * v[b])
        ii = np.where(v == -1)[0]
        # if there are fewer than 10 vertices, we abort; not enough to accurately
        # determine the values:
        if len(ii) == 0: return []
        (a,b) = (a[ii],b[ii])
        # Okay, get the average distance at which the crossing occurs
        d = np.mean([d0[a], d0[b]], axis=0)
        d = np.mean(d)
        # Before we find the vertices within this distance, we apply the masks so
        # that we exclude only those vertices not in these masks
        ii = np.setdiff1d(hem.indices, mask_ii)
        d0[ii] = np.nan # nans are never less than d
        ii = np.where(ny.util.nanle(d0, d))[0]
        return ii #if len(ii) >= 10 else []
    def generate_DROIs(self, sid, h, dist_prop, val_prop, ref_val, val_diff,
                      masks=None, distance_masks=None, inv_prop=None):
        '''
        generate_DROIs(sid, h, dist_prop, val_prop, ref_val, val_diff) yields a mask
          for the given subject/hemisphere that is based on the distance (which is
          measured by the dist_prop) required for the given val_prop to go from
          ref_val (the presumed value at distance = 0) to ref_val +/- val_diff.
        '''
        return _generate_DROIs(self.subjects[sid], h, dist_prop, ref_val, val_diff,
                               masks=masks, distance_masks=distance_masks, inv_prop=inv_prop)
    @staticmethod
    def _vertical_DROI_from_ventral_dorsal(vnt, drs):
        '''
        Given vnt and drs return values from the generate_DROI_data() function,
        yields the 'vertical' ROI (the V1 parts of vnt and drs) of the combined ROI.
        '''
        import neuropythy as ny, numpy as np, six
        (iiv, iid) = [np.where(q['visual_area'] == 1)[0] for q in (vnt,drs)]
        if len(iiv) == 0: return {k:v[iid] for (k,v) in six.iteritems(drs)}
        if len(iid) == 0: return {k:v[iiv] for (k,v) in six.iteritems(vnt)}
        res = {}
        for (k,v) in six.iteritems(vnt):
            res[k] = np.concatenate([v[iiv], drs[k][iid]])
        return res
    @staticmethod
    def _generate_subject_DROI_boundary_data(sub, h, paradigm, angle_delta,
                                             min_variance_explained=0,
                                             eccentricity_range=(0,7), surface_area='midgray'):
        import neuropythy as ny, numpy as np, copy, six
        CLS = VisualPerformanceFieldsDataset
        paradigm = paradigm.lower()
        erng = eccentricity_range
        minvexpl = min_variance_explained
        if paradigm == 'vertical':
            # This is the weird case: we handle it separately: just run the function
            # for both ventral and dorsal and concatenate the V1 parts
            (vnt,drs) = [CLS._generate_subject_DROI_boundary_data(sub, h, para, angle_delta,
                                                          min_variance_explained=minvexpl,
                                                          eccentricity_range=erng,
                                                          surface_area=surface_area)
                         for para in ['ventral', 'dorsal']]
            f = CLS._vertical_DROI_from_ventral_dorsal
            return f(vnt, drs)
        # Get the hemisphere:
        hem = sub.hemis[h]
        # Setup masks and handle eccentricity_range
        if erng in [None,Ellipsis]: erng = (0,7)
        if pimms.is_number(erng): erng = (0,erng)
        masks  = [('inf_eccentricity', erng[0], erng[1])]
        dmasks = [('prf_variance_explained', min_variance_explained, 1)]
        # Get the inferred angle (we'll need this)
        angle = np.abs(hem.prop('inf_polar_angle'))
        # We setup two different ROIs for V1/V2 or for D/V then we join them; this
        # is because we aren;t 100% confident that the boundaries are drawn in the
        # right place, but this should let the ROI grow appropriately on either side
        masks  = (copy.copy(masks),  masks)
        dmasks = (copy.copy(dmasks), dmasks)
        if paradigm == 'ventral':
            dprop = 'ventral_distance'
            xprop = 'dorsal_distance'
            ref_angle = 0
            for m in masks: m.append(('inf_polar_angle', -0.1, 90))
            masks[0].append(('inf_visual_area', 1))
            masks[1].append(('inf_visual_area', 2))
        elif paradigm == 'dorsal':
            dprop = 'dorsal_distance'
            xprop = 'ventral_distance'
            ref_angle = 180
            for m in masks: m.append(('inf_polar_angle', 90, 180.1))
            masks[0].append(('inf_visual_area', 1))
            masks[1].append(('inf_visual_area', 2))
        elif paradigm == 'horizontal':
            dprop = 'horizontal_distance'
            xprop = None
            ref_angle = 90
            for m in masks: m.append(('inf_visual_area', 1))
            masks[0].append(('inf_polar_angle', 90, 180.1))
            masks[1].append(('inf_polar_angle', -0.1,  90))
        else: raise ValueError('unrecognized paradigm: %s' % (paradigm,))
        # Get the indices
        ii = []
        for (m, dm) in zip(masks, dmasks):
            kk = CLS._generate_DROIs(sub, h, dprop, 'prf_polar_angle', ref_angle, angle_delta,
                                     masks=m, distance_masks=dm, inv_prop=xprop)
            ii = np.union1d(ii, kk)
        ii = np.array(ii, dtype='int')
        # Grab the other data:
        if pimms.is_str(surface_area) and not surface_area.endswith('_surface_area'):
            surface_area = surface_area + '_surface_area'
        surface_area = hem.property(surface_area)
        sa = surface_area[ii]
        th = hem.prop('thickness')[ii]
        vl = sa * th
        return {'surface_area_mm2': sa, 'mean_thickness_mm':th,
                'volume_mm3': vl,       'indices': ii,
                'visual_area': hem.prop('inf_visual_area')[ii]}
    def generate_subject_DROI_boundary_data(self, sid, h, paradigm, angle_delta,
                                            min_variance_explained=0,
                                            eccentricity_range=(0,7), surface_area='midgray'):
        '''
        generate_subject_DROI_boundary_data(sid, h, paradigm, delta) yields a dict of data about the
          given distance ROI; the ROI is defined by the paradigm, which must be one of 'vertical',
          'horizontal', 'ventral', or 'dorsal', and angle_delta, which is the distance in polar
          angle degrees from the given boundary.
        '''
        f = VisualPerformanceFieldsDataset._generate_subject_DROI_boundary_data
        return f(self.subjects[sid], h, paradigm, angle_delta, min_variance_explained=0,
                 eccentricity_range=(0,7), surface_area='midgray')
    @staticmethod
    def _generate_subject_DROI_data(sub, h, angle_delta, results='summary',
                                    eccentricity_range=(0,7), surface_area='midgray'):
        import neuropythy as ny, numpy as np, six
        results = results.lower()
        f = VisualPerformanceFieldsDataset._generate_subject_DROI_boundary_data
        (vnt,drs,hrz) = [f(sub, h, para, angle_delta,
                           eccentricity_range=eccentricity_range,
                           surface_area=surface_area)
                         for para in ['ventral','dorsal','horizontal']]
        # we don't need to run vertical because we can derive it from the other measures:
        ver = VisualPerformanceFieldsDataset._vertical_DROI_from_ventral_dorsal(vnt, drs)
        # depending on the results arg, we return these or their summaries
        res = {'vertical': ver, 'horizontal': hrz, 'ventral': vnt, 'dorsal': drs}
        if results == 'summary':
            fix = {'surface_area_mm2': np.nansum,
                   'volume_mm3': np.nansum,
                   'mean_thickness_mm': lambda x: np.nan if len(x) == 0 else np.nanmean(x)}
            return {k: {kk: fix[kk](vv) for (kk,vv) in six.iteritems(v) if kk in fix}
                    for (k,v) in six.iteritems(res)}
        else:
            return res
    def generate_subject_DROI_data(self, sid, h, angle_delta, results='summary',
                                   eccentricity_range=(0,7), surface_area='midgray'):
        '''
        generate_subject_DROI_data(sid, h, angle_delta) yields distance-based ROI data for
          the set of ROIs (ventral, dorsal, horizontal, vertical) for the given
          subject and hemisphere.
        '''
        return VisualPerformanceFieldsDataset._generate_subject_DROI_data(
            self.subjects[sid], h, angle_delta, results=results,
            eccentricity_range=eccentricity_range, surface_area=surface_area)
    @staticmethod
    def _generate_subject_DROI_table(subjects, sid, angles=None, eccens=None):
        import neuropythy as ny, numpy as np, six
        if angles in (None, Ellipsis, 'fine', 'all'):
            angles = VisualPerformanceFieldsDataset.roi_allangles
        elif angles == 'coarse':
            angles = VisualPerformanceFieldsDataset.roi_angles
        if eccens in (None, Ellipsis, 'all'):
            eccens = VisualPerformanceFieldsDataset.roi_eccens
        tbl = ny.auto_dict(None, [])
        sub = subjects[sid]
        f = VisualPerformanceFieldsDataset._generate_subject_DROI_data
        for h in ['lh','rh']:
            # go through the eccen ranges and angles:
            for (ang,erng) in [(a,e) for a in angles for e in eccens]:
                # get all the summary data:
                alldat = f(sub, h, ang, eccentricity_range=erng)
                for (para,dat) in six.iteritems(alldat):
                    # append to the appropriate columns:
                    tbl['sid'].append(sid)
                    tbl['hemisphere'].append(h)
                    tbl['boundary'].append(para)
                    tbl['min_eccentricity_deg'].append(erng[0])
                    tbl['max_eccentricity_deg'].append(erng[1])
                    tbl['angle_delta_deg'].append(ang)
                    for (k,v) in six.iteritems(dat):
                        tbl[k].append(v)
        return ny.to_dataframe(tbl)
    def generate_subject_DROI_table(self, sid, angles=None, eccens=None):
        """
        Calculate the distance-based ROIs for a single subject; indended for use with
        multiprocessing. This function will load the subject's data instead of running
        the calculation if the relevant data-file exists.
        """
        f = VisualPerformanceFieldsDataset._generate_subject_DROI_table
        return f(self.subjects, sid, angles=angles, eccens=eccens)
    @staticmethod
    def _generate_subject_DROI_details(subjects, sid, h, eccentricity_range=None):
        from neuropythy.util import curry
        import six, pyrsistent as pyr, numpy as np
        paradigm_order = ['dorsal', 'ventral', 'vertical', 'horizontal']
        roi_eccens = VisualPerformanceFieldsDataset.roi_eccens
        roi_angles = VisualPerformanceFieldsDataset.roi_angles
        e = eccentricity_range
        if e is None or e is Ellipsis:
            e = list(VisualPerformanceFieldsDataset.roi_eccens)
        if pimms.is_list(e) and all(pimms.is_tuple(q) for q in e):
            f = VisualPerformanceFieldsDataset._generate_subject_DROI_details
            res = [f(subjects, sid, h, eccentricity_range=q) for q in e]
            q = res[0]
            def _merge(p, a):
                r = {}
                for k in six.iterkeys(q[p][a]):
                    u = [u[p][a][k] for u in res if len(u[p][a][k]) > 0]
                    if len(u) == 0: u = np.asarray([], dtype=np.float)
                    else: u = np.concatenate(u)
                    r[k] = u
                return pyr.pmap(r)
            return pyr.pmap({k: pimms.lmap({a: curry(_merge, k, a) for a in roi_angles})
                             for k in paradigm_order})
        f0 = VisualPerformanceFieldsDataset._generate_subject_DROI_data
        f = lambda sid,h,k: f0(subjects[sid], h, k, eccentricity_range=e, results='all')
        lm0 = pimms.lmap({k: curry(f, sid, h, k) for k in roi_angles})
        pfn = lambda p: pimms.lmap({k:curry(lambda k:lm0[k][p], k) for k in roi_angles})
        return pimms.lmap({p: curry(pfn, p) for p in paradigm_order})
    def generate_subject_DROI_details(self, sid, h, eccentricity_range=None):
        """
        Calculate the distance-based ROI details for a single subject. Details contain
        similar data as the subject's DROI table, but is in a nested dictionary format.
        """
        f = VisualPerformanceFieldsDataset._generate_subject_DROI_details
        return f(self.subjects, sid, h, eccentricity_range=eccentricity_range)
    def generate_DROI_tables(self, nprocs=None, printstatus=False):
        '''
        generate_DROI_tables() recalculates the set of distance-based ROIs for each subject
        based on the data in the inferred maps and pRFs.
        '''
        import neuropythy as ny, os, six, multithreading as mp
        drois = {}
        subject_list = self.subject_list
        nsubs = len(subject_list)
        f = VisualPerformanceFieldsDataset._generate_subject_DROI_tables
        f = curry(f, self.subjects)
        if nprocs is None or nprocs is Ellipsis: nprocs = mp.cpu_count()
        if nprocs < 2: nprocs = 1
        if nprocs == 1:
            for ii in np.arange(0, nsubs, nprocs):
                mx = np.min([len(subject_list), ii + nprocs])
                if printstatus:
                    print("%2d - %3d%%" % ((ii + 1)/nsubs * 100, mx/nsubs * 100))
                pool = mp.Pool(nprocs)
                sids = subject_list[ii:ii+nprocs]
                tbls = pool.map(f, sids)
                pool.close()
                # add these data into the overall roi table
                for (sid,tbl) in zip(sids,tbls):
                    drois[sid] = tbl
        else:
            for (ii,sid) in enumerate(subject_list):
                if printstatus and ii % 10 == 9:
                    print('%3d (%d), %4.1f%%' % (ii, sid, ii/len(subject_list)*100.0))
                drois[sid] = f(sid)
        return pimms.persist(drois)
    def generate_DROI_table(self):
        '''
        generate_DROI_table() is a method that recalculates the DROI summary table from the
        individual subject DROIs.
        '''
        import pandas
        df = pandas.concat(list(self.DROI_tables.values()))
        df.set_index(['sid','hemisphere'])
        return df
    
    # Distance-based ROIs ##########################################################################
    # The distance-based/wedge ROIs themselves.
    @pimms.value
    def DROI_details(subjects):
        '''
        DROI_details is a nested-dictionary structure of the various DROI details of each subject
        and hemisphere.
        '''
        import neuropythy as ny, os, six
        from neuropythy.util import curry
        f = curry(VisualPerformanceFieldsDataset._generate_subject_DROI_details, subjects)
        m = {sid: pimms.lmap({h: curry(f, sid, h) for h in ['lh','rh']})
             for sid in six.iterkeys(subjects)}
        return pimms.persist(m)
    @pimms.value
    def DROI_tables(subject_list, pseudo_path):
        '''
        DROI_tables (distance-based regions of interest) is a dictionary of ROIS used in the visual
        performance field project. DROI_tables[sid] is a dataframe of the ROI-data for the subject
        with ID sid.
        '''
        import neuropythy as ny, os, six
        # Load one subject.
        def _load_DROI(sid):
            # get a subject-specific cache_path
            cpath = pseudo_path.local_path('DROIs', '%s.csv' % (sid,))
            return ny.load(cpath)
        return pimms.lmap({sid: ny.util.curry(_load_DROI, sid) for sid in subject_list})
    @pimms.value
    def DROI_table(pseudo_path):
        '''
        DROI_table (distance-based ROI table) is a dataframe summarizing all the data from all the
        hemispheres and all the distance-based wedge ROIs used in the visual performance fields
        project.
        '''
        import neuropythy as ny
        df = ny.load(pseudo_path.local_path('DROI_table.csv'))
        df.set_index(['sid','hemisphere'])
        return df
    @pimms.value
    def DROI_summary(DROI_table):
        '''
        DROI_summary (distance-based ROI summary) is a nested-dictionary data structure that
        provides easily-plottable summaries of the DROI_table. A value DROI_summary[b][k][u][s]
        corresponds to the boundary b ('ventral', 'dorsal', 'vertical', or 'horizontal'), the
        key k ('surface_area_mm2', 'mean_thickness_mm', or 'volume_mm3'), angle bin u (where 
        u = 0, 1, 2, 3, 4 indicates 0-10, 10-20, 20-30, 30-40, 40-50 degrees away from the relevant
        boundary), and subject s. Note that subject number is not necessarily consistent across
        boundaries and keys as some subjects have no vertices in small ROIs and thus get excluded
        from the relevant collections. All data are collapsed across eccentricities from 1-6
        degrees.
        '''
        import neuropythy as ny, numpy as np
        def _dfsel(df, k, ang, emns=[1,2,3,4,5], emxs=[2,3,4,5,6]):
            tbls = [ny.util.dataframe_select(df, angle_delta_deg=ang,
                                             min_eccentricity_deg=mn,
                                             max_eccentricity_deg=mx)
                    for (mn,mx) in zip(emns,emxs)]
            tbls = [tbl[['sid','hemisphere',k]] for tbl in tbls]
            tbl = tbls[0]
            for t in tbls[1:]:
                tt = tbl.merge(t, on=['sid','hemisphere'])
                tt[k] = tt[k+'_x'] + tt[k+'_y']
                tbl = tt[['sid','hemisphere',k]]
            tl = tbl.loc[tbl['hemisphere'] == 'lh']
            tr = tbl.loc[tbl['hemisphere'] == 'rh']
            tt = tl.merge(tr, on='sid')
            tt = tt.sort_values('sid')
            return tt[k+'_x'].values + tt[k+'_y'].values
        dat = {
            para: {
                k: tuple([x0[np.isfinite(x0)]
                          for ang in VisualPerformanceFieldsDataset.roi_angles
                          for x0 in [_dfsel(df, k, ang)]])
                for k in ['surface_area_mm2', 'mean_thickness_mm', 'volume_mm3']}
            for para in ['horizontal','vertical','dorsal','ventral']
            for df in [ny.util.dataframe_select(DROI_table, boundary=para)]}
        return pimms.persist(dat)
    
neuropythy.datasets.core.add_dataset('visual_performance_fields',
                                     lambda:VisualPerformanceFieldsDataset())
