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
neuropythy.datasets.core.add_dataset('visual_performance_fields',
                                     lambda:VisualPerformanceFieldsDataset())
