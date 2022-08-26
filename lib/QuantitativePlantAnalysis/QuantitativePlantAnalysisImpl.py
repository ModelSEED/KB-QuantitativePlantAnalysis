# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import sys
import os
import uuid
import json
from os.path import exists
import pandas as pd
from QuantitativePlantAnalysis.sdkhelper import SDKHelper
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace

def format_numbers(number):
    if number < 0.0001:
        return "{:.6f}".format(number)
    if number < 0.001:
        return "{:.5f}".format(number)
    if number < 0.01:
        return "{:.4f}".format(number)
    if number < 0.1:
        return "{:.3f}".format(number)
    return "{:.2f}".format(number)
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
        params = SDKHelper.validate_args(params,["workspace"],{
            "biomass_composition":[],
            "nitrogen_source":[],
            "hemicellulose_fraction":[],
            "monomer_lignin_fraction":[],
            "organic_acid_fraction":[],
            "transport_costs":[]
        })
        #Constants
        ATPyld,NADPHy,hexose,Nmass = 28.0,11.5,180.1559,14.0067
        NfracP,NH4upt,NO3upt = 0.1709,0.0,3. * hexose / (Nmass * ATPyld)
        NO3asm,N2asm = 4. * hexose / (Nmass * NADPHy),12. * hexose / (Nmass * ATPyld)
        N2tran,Xcell,Xstarc,Xsugar = 3. * hexose / (Nmass * ATPyld),1.111,1.149,1.000
        XhemiC,XhemiD,XhemiG,XlignH = 1.228,1.230,1.265,2.558
        XlignG,XlignS,Xlipid,Xprotn = 2.605,2.638,2.866,1.832
        XoaAC,XoaMO,XoaOA,Xminrl,Xsporopollenin,Xsuberin = 0.794,0.642,4.002,0.15,3.2,3.3
        result_table = pd.DataFrame({})
        
        if len(params["biomass_composition"]) == 0:
            params["biomass_composition"] = [{
                "cellulose":0.33,"sporopollenin":0,"suberin":0,"hemicellulose":0.30,"starch":0.30,"sugars":0.04,
                "lignin":0.05,"lipid":0.03,"protein":0.13,"organic_acid":0.03,
                "minerals":0.05
            }]
        for bc in params["biomass_composition"]:
            #Setting variables and normalizing
            Cellul = bc["cellulose"]
            Sporopollenin = bc["sporopollenin"]
            Suberin = bc["suberin"]
            Hemice = bc["hemicellulose"]
            Starch = bc["starch"]
            Sugars = bc["sugars"]
            Lignin = bc["lignin"]
            Lipids = bc["lipid"]
            Proten = bc["protein"]
            OrAcid = bc["organic_acid"]
            Minerl = bc["minerals"]
            Total = Cellul + Hemice + Starch + Sugars + Lignin + Lipids + Proten + OrAcid + Minerl + Sporopollenin + Suberin
            Cellul = Cellul / Total
            Hemice = Hemice / Total
            Starch = Starch / Total
            Sugars = Sugars / Total
            Lignin = Lignin / Total
            Lipids = Lipids / Total
            Proten = Proten / Total
            OrAcid = OrAcid / Total
            Minerl = Minerl / Total
            Sporopollenin = Sporopollenin / Total
            Suberin = Suberin / Total
            GRcell = Xcell  * Cellul
            GRstar = Xstarc * Starch
            GRsugr = Xsugar * Sugars
            GRlipd = Xlipid * Lipids
            GRprot = Xprotn * Proten
            GRminl = Xminrl * Minerl
            GRsporopollenin = Xsporopollenin * Sporopollenin
            GRsuberin = Xsuberin * Suberin
            
            NewTot = Cellul + Hemice + Starch + Sugars + Lignin + Lipids + Proten + OrAcid + Minerl + Sporopollenin + Suberin
            
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
                    params["hemicellulose_fraction"] = [{"hemic":0.6,"hemid":0.2,"hemig":0.3}]
                for hf in params["hemicellulose_fraction"]:
                    #Setting variables and normalizing
                    HemiC  = hf["hemic"]
                    HemiD  = hf["hemid"]
                    HemiG  = hf["hemig"]
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
                            GluReq = GRNacq + GRcell + GRstar + GRsporopollenin + GRsuberin + GRsugr + GRlipd + GRprot + GRminl + GRhemi + GRlign + GRoa
                            if len(params["transport_costs"]) == 0:
                                params["transport_costs"] = [{"cost_selection":"trans_phloem","fraction_starch":0.5,"membrane_crossings":1}]
                            for tc in params["transport_costs"]:
                                CostOptions = tc["cost_selection"]
                                FractionStarch = tc["fraction_starch"]
                                MembraneCrossings = tc["membrane_crossings"]
                                
                                ATPstarch = 1.0
                                ATPmembrane = 1.0
                                ATPyieldResp = 28.0
                                StaMobil = ATPstarch * FractionStarch * GluReq
                                Sucrose = 1.0 - FractionStarch
                                SucAmount = GluReq - StaMobil
                                TransCost = 0
                                if CostOptions == "trans_phloem":
                                    TransCost = ((MembraneCrossings * ATPmembrane * GluReq / 2.0) + StaMobil) / ATPyieldResp
                                elif CostOptions == "trans_no_phloem":
                                    TransCost = StaMobil / ATPyieldResp
                                
                                current_output = {}
                                current_output["Biomass composition"] = "<table>"\
                                    '<tr><th>Biomass constituent</th><th>Fraction</th><th>GluReq</th></tr>'+ \
                                    '<tr><td>Cellulose</td><td>'+format_numbers(Cellul)+"</td><td>"+format_numbers(GRcell)+"</td></tr>"+ \
                                    '<tr><td>Hemicelluloses</td><td>'+format_numbers(Hemice)+"</td><td>"+format_numbers(GRhemi)+"</td></tr>"+ \
                                    '<tr><td>Starch</td><td>'+format_numbers(Starch)+"</td><td>"+format_numbers(GRstar)+"</td></tr>"+ \
                                    '<tr><td>Sugars</td><td>'+format_numbers(Sugars)+"</td><td>"+format_numbers(GRsugr)+"</td></tr>"+ \
                                    '<tr><td>Lignins</td><td>'+format_numbers(Lignin)+"</td><td>"+format_numbers(GRlign)+"</td></tr>"+ \
                                    '<tr><td>Lipids</td><td>'+format_numbers(Lipids)+"</td><td>"+format_numbers(GRlipd)+"</td></tr>"+ \
                                    '<tr><td>Protein</td><td>'+format_numbers(Proten)+"</td><td>"+format_numbers(GRprot)+"</td></tr>"+ \
                                    '<tr><td>Organic Acids</td><td>'+format_numbers(OrAcid)+"</td><td>"+format_numbers(GRoa)+"</td></tr>"+ \
                                    '<tr><td>Minerals</td><td>'+format_numbers(Minerl)+"</td><td>"+format_numbers(GRminl)+"</td></tr>"+ \
                                    '<tr><td>N uptake/assimilation</td><td></td><td>'+format_numbers(GRNacq)+"</td></tr>"+ \
                                    '<tr><td>TOTAL</td><td>'+format_numbers(NewTot)+"</td><td>"+format_numbers(GluReq)+"</td></tr></table>"
                                current_output["Nitrogen source<br>(g/g plant N)"] = "<table>"\
                                    '<tr><td>Plant N [estimated]</td><td>'+format_numbers(PlantN)+"</td></tr>"+ \
                                    '<tr><td>N from NH4-N</td><td>'+format_numbers(N_NH4)+"</td></tr>"+ \
                                    '<tr><td>N from NO3-N</td><td>'+format_numbers(N_NO3)+"</td></tr>"+ \
                                    '<tr><td>N from  N2-N</td><td>'+format_numbers(N_N2)+"</td></tr>"+ \
                                    '<tr><td>N assimilation cost</td><td>'+format_numbers(GRNacq)+"</td></tr></table>"
                                current_output["Hemicellulose fraction<br>(g/g plant)"] = "<table>"\
                                    '<tr><td>Hemi conifer</td><td>'+format_numbers(HemiC)+"</td></tr>"+ \
                                    '<tr><td>Hemi dicot</td><td>'+format_numbers(HemiD)+"</td></tr>"+ \
                                    '<tr><td>Hemi grass</td><td>'+format_numbers(HemiG)+"</td></tr></table>"
                                current_output["Monomer lignin fraction<br>(g/g plant)"] = "<table>"\
                                    '<tr><td>Coumaryl</td><td>'+format_numbers(Coumrl)+"</td></tr>"+ \
                                    '<tr><td>Coniferyl</td><td>'+format_numbers(Conifr)+"</td></tr>"+ \
                                    '<tr><td>Sinapyl</td><td>'+format_numbers(Sinapl)+"</td></tr></table>"
                                current_output["Organic acid fraction<br>(g/g plant)"] = "<table>"\
                                    '<tr><td>Aconitic citric</td><td>'+format_numbers(AcoCit)+"</td></tr>"+ \
                                    '<tr><td>Malic oxaloacetic</td><td>'+format_numbers(MalOxa)+"</td></tr>"+ \
                                    '<tr><td>Oxalic</td><td>'+format_numbers(Oxalic)+"</td></tr></table>"
                                current_output["Transport costs"] = "<table>"\
                                    '<tr><td>Cost options</td><td>'+CostOptions+"</td></tr>"+ \
                                    '<tr><td>Fraction starch</td><td>'+format_numbers(FractionStarch)+"</td></tr>"+ \
                                    '<tr><td>Membrane crossings</td><td>'+str(MembraneCrossings)+"</td></tr></table>"
                                current_output["Growth yield results"] = "<table>"\
                                    '<tr><td>Glucose req<br>(g/g plant)</td><td>'+format_numbers(GluReq)+"</td></tr>"+ \
                                    '<tr><td>Transport req<br>(g/g plant)</td><td>'+format_numbers(TransCost)+"</td></tr>"+ \
                                    '<tr><td>Total req<br>(g/g plant)</td><td>'+format_numbers(TransCost+GluReq)+"</td></tr>"+ \
                                    '<tr><td>Yield req<br>(g plant/g glucose)</td><td>'+format_numbers(1./(TransCost+GluReq))+"</td></tr></table>"
                                result_table = result_table.append(current_output, ignore_index = True)
        column_list = ["Biomass composition","Nitrogen source<br>(g/g plant)","Hemicellulose fraction<br>(g/g plant)","Monomer lignin fraction<br>(g/g plant)","Organic acid fraction<br>(g/g plant)","Growth yield results"]
        html_data = f"""
            <html>
            <header>
                <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
            </header>
            <body>
            {result_table.to_html(columns=column_list,escape=False,notebook=False,table_id="table",index=False,justify="left")}
            <script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
            <script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
            <script>
                $(document).ready( function () {{
                    $('#table').DataTable({{
                        // paging: false,    
                        // scrollY: 400,
                    }});
                }});
            </script>
            </body>
            </html>
            """
        report_name = str(uuid.uuid4())
        workspace = None
        if isinstance(params["workspace"], str):
            workspace = params["workspace"]
        else:
            workspace = str(params["workspace"])
        html_report_folder = os.path.join(self.shared_folder, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        with open(os.path.join(html_report_folder, 'index.html'), 'w') as f:
            f.write(html_data)
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
        repout = report.create_extended_report(report_params)
        output = {"report_name":report_name,"report_ref":repout["ref"],'workspace_name':workspace}
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
