# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import sys
import os
import uuid
from os.path import exists
import pandas as pd
from QuantitativePlantAnalysis.sdkhelper import SDKHelper
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace
#END_HEADER


class QuantitativePlantAnalysis:
    '''
    Module Name:
    QuantitativePlantAnalysis

    Module Description:
    A KBase module: QuantitativePlantAnalysis
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def compute_plant_biomass_yield(self, ctx, params):
        """
        Compute plant biomass yield from glucose
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN compute_plant_biomass_yield
        #Processing parameters
        SDKHelper.validate_args(params,["workspace"],{
            "biomass_composition":[],
            "nitrogen_source":[],
            "hemicellulose_fraction":[],
            "monomer_lignin_fraction":[],
            "organic_acid_fraction":[]
        })
        #Constants
        ATPyld,NADPHy,hexose,Nmass = 28.0,11.5,180.1559,14.0067
        NfracP,NH4upt,NO3upt = 0.1709,0.0,3. * hexose / (Nmass * ATPyld)
        NO3asm,N2asm = 4. * hexose / (Nmass * NADPHy),12. * hexose / (Nmass * ATPyld)
        N2tran,Xcell,Xstarc,Xsugar = 3. * hexose / (Nmass * ATPyld),1.111,1.149,1.000
        XhemiC,XhemiD,XhemiG,XlignH = 1.228,1.230,1.265,2.558
        XlignG,XlignS,Xlipid,Xprotn = 2.605,2.638,2.866,1.832
        XoaAC,XoaMO,XoaOA,Xminrl = 0.794,0.642,4.002,0.15
        #"Monomer lignin fraction":None,"Organic acid fraction":None,"Hemicellulose fraction":None
        default_output = {"Biomass composition":None,"Nitrogen source":None,
                          "Growth yield results":None}
        result_table = pd.DataFrame({})
        
        if len(params["biomass_composition"]) == 0:
            params["biomass_composition"] = [{
                "cellulose":0.33,"hemicellulose":0.30,"starch":0.30,"sugars":0.04,
                "lignin":0.05,"lipid":0.03,"protein":0.13,"organic_acid":0.03,
                "minerals":0.05
            }]
        for bc in params["biomass_composition"]:
            #Setting variables and normalizing
            Cellul = bc["cellulose"]
            Hemice = bc["hemicellulose"]
            Starch = bc["starch"]
            Sugars = bc["sugars"]
            Lignin = bc["lignin"]
            Lipids = bc["lipid"]
            Proten = bc["protein"]
            OrAcid = bc["organic_acid"]
            Minerl = bc["minerals"]
            Total = Cellul + Hemice + Starch + Sugars + Lignin + Lipids + Proten + OrAcid + Minerl
            Cellul = Cellul / Total
            Hemice = Hemice / Total
            Starch = Starch / Total
            Sugars = Sugars / Total
            Lignin = Lignin / Total
            Lipids = Lipids / Total
            Proten = Proten / Total
            OrAcid = OrAcid / Total
            Minerl = Minerl / Total
            GRcell = Xcell  * Cellul
            GRstar = Xstarc * Starch
            GRsugr = Xsugar * Sugars
            GRlipd = Xlipid * Lipids
            GRprot = Xprotn * Proten
            GRminl = Xminrl * Minerl
            
            NewTot = Cellul + Hemice + Starch + Sugars + Lignin + Lipids + Proten + OrAcid + Minerl
            
            if len(params["nitrogen_source"]) == 0:
                params["nitrogen_source"] = [{"nh4":0.3,"no3":0.69,"n2":0.01}]
            for ns in params["nitrogen_source"]:
                #Setting variables and normalizing
                NH4 = ns["nh4"]
                NO3 = ns["no3"]
                N2 = ns["n2"]
                Total = NH4 + NO3 + N2
                NH4    = NH4 / Total
                NO3    = NO3 / Total
                N2     = N2  / Total
                
                PlantN = NfracP * Proten
                N_NH4  = NH4    * PlantN
                N_NO3  = NO3    * PlantN
                N_N2   = N2     * PlantN
                GRNacq = NH4upt * N_NH4 + (NO3upt + NO3asm) * N_NO3 + (N2asm + N2tran) * N_N2
                
                if len(params["hemicellulose_fraction"]) == 0:
                    params["hemicellulose_fraction"] = [{"HemiC":0.6,"HemiD":0.2,"HemiG":0.3}]
                for hf in params["hemicellulose_fraction"]:
                    #Setting variables and normalizing
                    HemiC  = hf["HemiC"]
                    HemiD  = hf["HemiD"]
                    HemiG  = hf["HemiG"]
                    Total = HemiC  + HemiD + HemiG
                    HemiC  = HemiC  / Total
                    HemiD  = HemiD  / Total
                    HemiG  = HemiG  / Total
                    
                    GRhemi = XhemiC * HemiC + XhemiD * HemiD + XhemiG * HemiD
                    GRhemi = GRhemi * Hemice
        
                    if len(params["monomer_lignin_fraction"]) == 0:
                        params["monomer_lignin_fraction"] = [{"coumaryl":0.4,"coniferyl":0.4,"sinapyl":0.2}]
                    for mlf in params["monomer_lignin_fraction"]:
                        #Setting variables and normalizing
                        Coumrl = mlf["coumaryl"]
                        Conifr = mlf["coniferyl"]
                        Sinapl = mlf["sinapyl"]
                        Total = Coumrl + Conifr + Sinapl
                        Coumrl = Coumrl / Total
                        Conifr = Conifr / Total
                        Sinapl = Sinapl / Total
                        
                        GRlign = XlignH * Coumrl + XlignG * Conifr + XlignS * Sinapl
                        GRlign = GRlign * Lignin
                        
                        if len(params["organic_acid_fraction"]) == 0:
                            params["organic_acid_fraction"] = [{"aconitic_citric":0.475,"malic_oxaloacetic":0.475,"oxalic":0.05}]
                        for oaf in params["organic_acid_fraction"]:
                            #Setting variables and normalizing
                            AcoCit = oaf["aconitic_citric"]
                            MalOxa = oaf["malic_oxaloacetic"]
                            Oxalic = oaf["oxalic"]
                            Total = MalOxa + AcoCit + Oxalic
                            AcoCit = AcoCit / Total
                            MalOxa = MalOxa / Total
                            Oxalic = Oxalic / Total
                            
                            GRoa   = XoaAC * AcoCit + XoaMO * MalOxa + XoaOA * Oxalic
                            GRoa   = GRoa * OrAcid
                            GluReq = GRNacq + GRcell + GRstar + GRsugr + GRlipd + GRprot + GRminl + GRhemi + GRlign + GRoa
                            
                            current_output = default_output.copy()
                            
                            current_output["Biomass composition"] = "<table>"\
                                '<tr><th>Biomass constituent</th><th>Fraction</th><th>GluReq</th></tr>'+ \
                                '<tr><td>Cellulose</td><td>'+str(Cellul)+"</td><td>"+str(GRcell)+"</td></tr>"+ \
                                '<tr><td>Hemicelluloses</td><td>'+str(Hemice)+"</td><td>"+str(GRhemi)+"</td></tr>"+ \
                                '<tr><td>Starch</td><td>'+str(Starch)+"</td><td>"+str(GRstar)+"</td></tr>"+ \
                                '<tr><td>Sugars</td><td>'+str(Sugars)+"</td><td>"+str(GRsugr)+"</td></tr>"+ \
                                '<tr><td>Lignins</td><td>'+str(Lignin)+"</td><td>"+str(GRlign)+"</td></tr>"+ \
                                '<tr><td>Lipids</td><td>'+str(Lipids)+"</td><td>"+str(GRlipd)+"</td></tr>"+ \
                                '<tr><td>Protein</td><td>'+str(Proten)+"</td><td>"+str(GRprot)+"</td></tr>"+ \
                                '<tr><td>Organic Acids</td><td>'+str(OrAcid)+"</td><td>"+str(GRoa)+"</td></tr>"+ \
                                '<tr><td>Minerals</td><td>'+str(Minerl)+"</td><td>"+str(GRminl)+"</td></tr>"+ \
                                '<tr><td>N uptake/assimilation</td><td>'+str(0.0)+"</td><td>"+str(GRNacq)+"</td></tr>"+ \
                                '<tr><td>TOTAL</td><td>'+str(NewTot)+"</td><td>"+str(GluReq)+"</td></tr></table>"
                            current_output["Nitrogen source"] = "<table>"\
                                '<tr><td>Plant N [estimated]</td><td>'+str(PlantN)+" (g/g plant)</td></tr>"+ \
                                '<tr><td>N from NH4-N</td><td>'+str(N_NH4)+" (g/g plant)</td></tr>"+ \
                                '<tr><td>N from NO3-N</td><td>'+str(N_NO3)+" (g/g plant)</td></tr>"+ \
                                '<tr><td>N from  N2-N</td><td>'+str(N_N2)+" (g/g plant)</td></tr>"+ \
                                '<tr><td>N assimilation cost</td><td>'+str(GRNacq)+" (g glu/g plant)</td></tr></table>"
                            current_output["Growth yield results"] = \
                                str(GluReq)+" (g glucose/g plant)<br>"+ \
                                str(1./GluReq)+" (g plant/g glucose)"
                            result_table = result_table.append(current_output, ignore_index = True)
        report_name = str(uuid.uuid4())
        workspace = None
        if isinstance(params["workspace"], str):
            workspace = params["workspace"]
        else:
            workspace = str(params["workspace"])
        output = {"report_name":report_name,"report_ref":workspace+"/"+report_name,'workspace_name':workspace,'ws':workspace}
        html_report_folder = os.path.join(self.shared_folder, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        with open(os.path.join(html_report_folder, 'index.html'), 'w') as f:
            f.write(result_table.to_html(escape=False,notebook=False))
        report_shock_id = self.dfu.file_to_shock({'file_path': html_report_folder,'pack': 'zip'})['shock_id']
        html_output = {
            'name' : 'index.html',
            'shock_id': report_shock_id
        }
        report_params = {
            'objects_created': [],
            'workspace_name': params["workspace"],
            'html_links': [{
                'name' : 'index.html',
                'shock_id': report_shock_id
            }],
            'direct_html_link_index': 0,
            'html_window_height': 700,
            'report_object_name': report_name
        }
        report = KBaseReport(self.callback_url, token=ctx['token'])
        report.create_extended_report(report_params)
        #END compute_plant_biomass_yield

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method compute_plant_biomass_yield return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
