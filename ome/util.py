from ome import base, settings
from ome.data import *
from ome.components import *

from IPython.display import HTML
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram

from IPython.display import display
from IPython.html import widgets as W
from IPython.utils import traitlets as T

import pandas as pd
import numpy as np
import cobra
import cPickle as pickle


model = pickle.load(open(settings.data_directory+'/models/iJO1366.pickle', "rb"))
#model = cobra.io.load_matlab_model(settings.data_directory+'/models/iJO1366')



def get_chip_peak_gene_expression(gene_name, factor, condition):

    #if gene_name = 'nan': return None

    chip_regulation = ome.query(ChipPeakGene).filter_by(target=factor, condition=condition, gene_name=str(gene_name)).all()

    if not chip_regulation: return None
    aa = AllAnalysis
    target = 'delta-'+factor[0].lower()+factor[1:]
    array_regulation = ome.query(aa).filter(or_(and_(aa.target1 == target, aa.target2 == 'wt'),\
                                                and_(aa.target2 == target, aa.target1 == 'wt')),\
                                            and_(aa.gene_name == gene_name, aa.condition1 == condition, aa.fdr < .05)).all()

    if len(array_regulation) > 1: return 'error'
    elif len(array_regulation) == 0: return None

    if array_regulation[0].target1 == 'wt': return {'fold_change':float(array_regulation[0].fold_change)}
    else: return {'fold_change': -1*float(array_regulation[0].fold_change)}


#def get_reaction_data(cobra_rxn_id, type=GeneExpressionData, **kwargs):



def get_regulation_data(cobra_rxn_id, verbose=True, as_dataframe=True):
    """This function takes a cobra_rxn_id and returns a json regulation structure, by default it works for chip peak gene expression
       and the numeric value always corresponds to positive(activation) and negative(repression)
    {'cobra_rxn_id': [{'value':1., 'type':'transcriptional', 'regulator':'component_1', 'gene':'gene_id_1', 'dataset':dataset_id_1},
                      {'value':-5, 'type':'transcriptional', 'regulator':'component_2', 'gene':'gene_id_2', 'dataset':dataset_id_1},
                      {'value':3., 'type':'translational', 'regulator':'component_id_3', 'mRNA':'rna_id_36', 'dataset':dataset_id_2},
                       ...
                      {'value':n, 'type':'allosteric', 'regulator':'component_id_n', 'protein':'protein_id_n', 'dataset':dataset_id_n}]}
    """
    try:
        cobra_rxn_id.id
        rxn = cobra_rxn_id
    except:
        rxn = model.reactions.get_by_id(cobra_rxn_id)

    cpge = ChIPPeakGeneExpression

    gpr_genes = [ome.query(Gene.name).filter_by(locus_id=gene.id).scalar() for gene in rxn.genes]

    gpr_regulation = ome.query(cpge.value, cpge.strain1, cpge.target, cpge.gene_name,cpge.diff_exp_id,
                               cpge.carbon_source, cpge.nitrogen_source, cpge.electron_acceptor).\
                                 filter(cpge.gene_name.in_(gpr_genes)).\
                                 group_by(cpge.value, cpge.strain1, cpge.target, cpge.gene_name,
                                          cpge.diff_exp_id, cpge.carbon_source, cpge.nitrogen_source,
                                          cpge.electron_acceptor).all()


    if verbose:
        data = [{'value':x.value, 'strain':x.strain1, 'type':'transcriptional', 'regulator':x.target, 'gene':x.gene_name,
                 'dataset':x.diff_exp_id, 'carbon_source':x.carbon_source, 'nitrogen_source':x.nitrogen_source,
                 'electron_acceptor':x.electron_acceptor} for x in gpr_regulation]
    else:
        data = [{'value':x.value, 'type':'transcriptional', 'regulator':x.target, 'gene':x.gene_name,
                 'dataset':x.diff_exp_id} for x in gpr_regulation]

    if as_dataframe:
    	try: return pd.DataFrame(data, columns=['dataset','regulator','strain','carbon_source',
    											'nitrogen_source','electron_acceptor','gene',
    											'type','value']).sort(['gene','strain'])
    	except: return []
    else: return data


def get_indirect_regulation_frame(rxn, factors=['ArcA','Fnr'], condition=None):
    multi_index = [[],[]]
    for factor in factors:
        for gene in rxn._genes.keys():
            if gene == 's0001': continue
            multi_index[0].append(factor)
            multi_index[1].append(ome.query(Gene.name).filter_by(bnum=gene.id).scalar())

    df = DataFrame(index=multi_index, columns=['value'])
    for factor,gene_name in df.index:
        aa = AllAnalysis
        target = 'delta-'+factor[0].lower()+factor[1:]
        array_regulation = ome.query(aa).filter(or_(and_(aa.target1 == target, aa.target2 == 'wt'),\
                                                    and_(aa.target2 == target, aa.target1 == 'wt')),\
                                                    and_(aa.gene_name == gene_name, aa.condition1 == condition, aa.fdr < .05)).all()

        if len(array_regulation) > 1: return 'error'
        elif len(array_regulation) == 0: continue

        if array_regulation[0].target1 == 'wt':
            df['value'].ix[factor,gene_name] = float(array_regulation[0].fold_change)
        else:
            df['value'].ix[factor,gene_name] = -1*float(array_regulation[0].fold_change)

    return df



def add_gene_group(name, genes, genome):
    session = Session()

    gene_group = session.get_or_create(GeneGroup, name = name)

    for gene in genes:
        if isinstance(gene, basestring):
            gene = session.query(Gene).filter(and_(Gene.genome_id == genome.id,
                                                   or_(Gene.name == gene,
                                                       Gene.locus_id == gene))).first()

        if gene: session.get_or_create(GeneGrouping, group_id=gene_group.id, gene_id=gene.id)

    session.flush()
    session.commit()
    session.close()



def write_genome_data_gff(genome_data_set_ids, function='avg'):
    session = base.Session()

    data_sets = session.query(DataSet).filter(DataSet.id.in_(genome_data_set_ids)).all()

    if function:
        vals = data_sets[0].name.split('_')
        #name = '_'.join(vals[0:5]+[6:]+[function]

    genbank_fasta_string = 'gi|'+genome.genbank_id+'|ref|'+genome.ncbi_id+'|'

    with open(settings.data_directory+'/annotation/'+genome.ncbi_id+'.gff', 'wb') as gff_file:

        for gene in session.query(components.Gene).all():

            info_string = 'gene_id "%s"; transcript_id "%s"; gene_name "%s";' % (gene.locus_id, gene.locus_id, gene.name)

            gff_string = '%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s\n' % (genbank_fasta_string, 'ome_db', 'exon', gene.leftpos,
                                                                                                         gene.rightpos,
                                                                                                         gene.strand,
                                                                                                         info_string)
            gff_file.write(gff_string)

    session.close()




def gene_heatmap(gene_list, analysis_type=GeneExpressionData, dataset_type='%',
                                                              strain1=ome.query(Strain).all(),
                                                              strain2=ome.query(Strain).all(),
                                                              environments1=ome.query(InVivoEnvironment).all(),
                                                              environments2=ome.query(InVivoEnvironment).all()):


    dataset_id_map = {}
    if analysis_type == GeneExpressionData:
        all_data = ome.query(analysis_type).filter(and_(analysis_type.gene_name.in_([g.name for g in gene_list]),
                                                        analysis_type.strain_id.in_([x.id for x in strain1]),
                                                        analysis_type.environment_id.in_([x.id for x in environments1]),
                                                        analysis_type.dataset_type.ilike(dataset_type))).all()

        datasets = set(['_'.join([x.dataset_type,x.strain,x.carbon_source,x.nitrogen_source,x.electron_acceptor]) for x in all_data])


        genes_data = pd.DataFrame(index=[g.name for g in gene_list], columns=datasets)

        for x in all_data:
            if dataset_type == 'array_experiment':
                genes_data.ix[x.gene_name]['_'.join([x.dataset_type,x.strain,x.carbon_source,x.nitrogen_source,x.electron_acceptor])] = x.value
            else:
                genes_data.ix[x.gene_name]['_'.join([x.dataset_type,x.strain,x.carbon_source,x.nitrogen_source,x.electron_acceptor])] = np.log(x.value)


    elif analysis_type == DifferentialGeneExpressionData:
        all_data = ome.query(analysis_type).\
                              filter(and_(analysis_type.gene_id.in_([g.id for g in gene_list]),
                                          analysis_type.strain1.in_([x.name for x in strain1]),
                                          analysis_type.strain2.in_([x.name for x in strain2]),
                                          analysis_type.environment_id_1.in_([x.id for x in environments1]),
                                          analysis_type.environment_id_2.in_([x.id for x in environments2]))).all()

        datasets = set([x.diff_name() for x in all_data])
        genes_data = pd.DataFrame(index=[g.name for g in gene_list], columns=datasets)

        for x in all_data:
            genes_data.ix[x.gene_name][x.diff_name()] = x.value


    elif analysis_type == ChIPPeakGeneExpression:
        all_data = ome.query(analysis_type).\
                              filter(and_(analysis_type.gene_id.in_([g.id for g in gene_list]),
                                          analysis_type.strain1.in_([x.name for x in strain1]),
                                          analysis_type.strain2.in_([x.name for x in strain2]),
                                          analysis_type.environment_id.in_([x.id for x in environments2]))).all()

        datasets = set(['_'.join([x.target,x.strain1+'/'+x.strain2,x.carbon_source,x.nitrogen_source,x.electron_acceptor]) for x in all_data])
        genes_data = pd.DataFrame(index=[g.name for g in gene_list], columns=datasets)

        for x in all_data:
            genes_data.ix[x.gene_name]['_'.join([x.target,x.strain1+'/'+x.strain2,x.carbon_source,x.nitrogen_source,x.electron_acceptor])] = x.value



    genes_data = genes_data.dropna(how='all').fillna(0.)
    genes_data = genes_data.replace([np.inf], 10.)
    genes_data = genes_data.replace([-np.inf], -10.)
    col_labels = list(genes_data.index)
    row_labels = list(datasets)


    heatmap_data = []
    for i,g in enumerate(genes_data.index):
        for j,c in enumerate(genes_data.columns):
            heatmap_data.append({"row": j+1, "col": i+1, "value": genes_data.ix[g][c]})


    dm = genes_data
    D1 = squareform(pdist(dm, metric='euclidean'))
    D2 = squareform(pdist(dm.T, metric='euclidean'))

    Y = linkage(D1, method='single')
    Z1 = dendrogram(Y, labels=dm.index)

    Y = linkage(D2, method='single')
    Z2 = dendrogram(Y, labels=dm.columns)

    hccol = Z1['leaves']
    hcrow = Z2['leaves']

    hccol = [x+1 for x in hccol]
    hcrow = [x+1 for x in hcrow]


    html_style = """
     <style>
      /* disable text selection */
      svg *::selection {
         background : transparent;
      }

      svg *::-moz-selection {
         background:transparent;
      }

      svg *::-webkit-selection {
         background:transparent;
      }
      rect.selection {
        stroke          : #333;
        stroke-dasharray: 4px;
        stroke-opacity  : 0.5;
        fill            : transparent;
      }

      rect.cell-border {
        stroke: #eee;
        stroke-width:0.3px;
      }

      rect.cell-selected {
        stroke: rgb(51,102,153);
        stroke-width:0.5px;
      }

      rect.cell-hover {
        stroke: #F00;
        stroke-width:0.3px;
      }

      text.mono {
        font-size: 10pt;
        font-family: Consolas, courier;
        fill: #aaa;
      }

      text.text-selected {
        fill: #000;
      }

      text.text-highlight {
        fill: #c00;
      }
      text.text-hover {
        fill: #00C;
      }
      #tooltip {
        position: fixed;
        width: 200px;
        height: auto;
        padding: 10px;
        background-color: white;
        -webkit-border-radius: 10px;
        -moz-border-radius: 10px;
        border-radius: 10px;
        -webkit-box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
        -moz-box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
        box-shadow: 4px 4px 10px rgba(0, 0, 0, 0.4);
        pointer-events: none;
      }

      #tooltip.hidden {
        display: none;
      }

      #tooltip p {
        margin: 0;
        font-family: sans-serif;
        font-size: 12px;
        line-height: 20px;
      }
      </style>


      <div id="tooltip" class="hidden">
          <p><span id="value"></p>
      </div>


      <select id="order">
      <option value="hclust">by cluster</option>
      <option value="probecontrast">by probe name and contrast name</option>
      <option value="probe">by probe name</option>
      <option value="contrast">by contrast name</option>
      <option value="custom">by log2 ratio</option>
      </select>
      </select>
      <div id="chart" style='overflow:auto; width:960px; height:auto;'></div>
      """



    #return hcrow,hccol,row_labels,col_labels,heatmap_data,dm
    return {'hcrow': hcrow, 'hccol': hccol, 'row_labels':row_labels,
                                            'col_labels':col_labels,
                                            'heatmap_data':heatmap_data,
                                            'maxval' : max([x['value'] for x in heatmap_data]),
                                            'minval' : min([x['value'] for x in heatmap_data]),
                                            'html_style': html_style}


class HeatmapWidget(W.DOMWidget):
    _view_name = T.Unicode('HeatmapView', sync=True)
    heatmap_data = T.List(sync=True)
    row_labels = T.List(sync=True)
    col_labels = T.List(sync=True)
    hcrow = T.List(sync=True)
    hccol = T.List(sync=True)
    minval = T.Float(sync=True)
    maxval = T.Float(sync=True)
    html_style = T.Unicode(sync=True)


heatmap_js = """
require.config({paths: {d3: "https://mpld3.github.io/js/d3.v3.min"}});

require(["widgets/js/widget", "d3"], function(WidgetManager, d3){
    var HeatmapView = IPython.DOMWidgetView.extend({

        render: function(){

            this.$el.append(this.model.get("html_style"));

            this.update();
        },
        update: function(){
            var margin = { top: 175, right: 10, bottom: 150, left: 400 },
                cellSize=18,
                 // - margin.top - margin.bottom,
                  //gridSize = Math.floor(width / 24),
                legendElementWidth = cellSize*2.5,
                colorBuckets = 21;
                //colors = ['#005824','#1A693B','#347B53','#4F8D6B','#699F83','#83B09B','#9EC2B3','#B8D4CB','#D2E6E3','#EDF8FB','#FFFFFF','#F1EEF6','#E6D3E1','#DBB9CD','#D19EB9','#C684A4','#BB6990','#B14F7C','#A63467','#9B1A53','#91003F'];
            var colors = ["#081d58","#162876","#253494","#23499E","#2253A3",#225ea8","#1F77B4","#1d91c0","#2FA3C2","#38ACC3","#41b6c4","#60C1BF","#7fcdbb","#91D4B9",#A3DBB7","#c7e9b4","#DAF0B2","#E3F4B1",#edf8b1","#F6FBC5","#ffffd9"];
            var rowLabel = this.model.get("row_labels");
            var colLabel = this.model.get("col_labels");
            var data = this.model.get("heatmap_data");
            var hcrow = this.model.get("hcrow");
            var hccol = this.model.get("hccol");

            var col_number = colLabel.length,
                row_number = rowLabel.length;

            var width = cellSize*col_number, // - margin.left - margin.right,
                height = cellSize*row_number;

            var maxval = JSON.parse(this.model.get("maxval"));
            var minval = JSON.parse(this.model.get("minval"));

            var colorScale = d3.scale.quantile()
                               .domain([ minval , 0, maxval])
                               .range(colors);

            this.svg = d3.select(this.el).append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            var svg = this.svg;

            var rowSortOrder=false;
            var colSortOrder=false;
            var rowLabels = svg.append("g")
                  .selectAll(".rowLabelg")
                  .data(rowLabel)
                  .enter()
                  .append("text")
                  .text(function (d) { return d; })
                  .attr("x", 0)
                  .attr("y", function (d, i) { return hcrow.indexOf(i+1) * cellSize; })
                  .style("text-anchor", "end")
                  .attr("transform", "translate(-25," + cellSize / 1.5 + ")")
                  .attr("class", function (d,i) { return "rowLabel mono r"+i;} )
                  .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
                  .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
                  .on("click", function(d,i) {rowSortOrder=!rowSortOrder; sortbylabel("r",i,rowSortOrder); this.$('#order').property("selectedIndex", 4).node().focus();;});

            var colLabels = svg.append("g")
                  .selectAll(".colLabelg")
                  .data(colLabel)
                  .enter()
                  .append("text")
                  .text(function (d) { return d; })
                  .attr("x", 0)
                  .attr("y", function (d, i) { return hccol.indexOf(i+1) * cellSize; })
                  .style("text-anchor", "left")
                  .attr("transform", "translate("+cellSize/2 + ",-25) rotate (-90)")
                  .attr("class",  function (d,i) { return "colLabel mono c"+i;} )
                  .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
                  .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
                  .on("click", function(d,i) {colSortOrder=!colSortOrder;  sortbylabel("c",i,colSortOrder); this.$('#order').property("selectedIndex", 4).node().focus();;});

            var heatMap = svg.append("g").attr("class","g3")
                  .selectAll(".cellg")
                  .data(data,function(d){return d.row+":"+d.col;})
                  .enter()
                  .append("rect")
                  .attr("x", function(d) { return hccol.indexOf(d.col) * cellSize; })
        .attr("y", function(d) { return hcrow.indexOf(d.row) * cellSize; })
        .attr("class", function(d){return "cell cell-border cr"+(d.row-1)+" cc"+(d.col-1);})
        .attr("width", cellSize)
        .attr("height", cellSize)
        .style("fill", function(d) { return colorScale(d.value); })
        /* .on("click", function(d) {
               var rowtext=d3.select(".r"+(d.row-1));
               if(rowtext.classed("text-selected")==false){
                   rowtext.classed("text-selected",true);
               }else{
                   rowtext.classed("text-selected",false);
               }
        })*/
        .on("mouseover", function(d){
               //highlight text
               d3.select(this).classed("cell-hover",true);
               d3.selectAll(".rowLabel").classed("text-highlight",function(r,ri){ return ri==(d.row-1);});
               d3.selectAll(".colLabel").classed("text-highlight",function(c,ci){ return ci==(d.col-1);});

               //Update the tooltip position and value
               d3.select("#tooltip")
                 .style("left", (d3.event.pageX+10) + "px")
                 .style("top", (d3.event.pageY-10) + "px")
                 .select("#value")
                 .text(rowLabel[d.row-1]+","+colLabel[d.col-1]+"\ndata:"+d.value);
               //Show the tooltip
               d3.select("#tooltip").classed("hidden", false);
        })
        .on("mouseout", function(){
               d3.select(this).classed("cell-hover",false);
               d3.selectAll(".rowLabel").classed("text-highlight",false);
               d3.selectAll(".colLabel").classed("text-highlight",false);
               d3.select("#tooltip").classed("hidden", true);
        })
        ;

  var legend = svg.selectAll(".legend")
      .data(d3.range(Math.floor(minval), Math.ceil(maxval)))
      .enter().append("g")
      .attr("class", "legend");

  var color_factor = Math.ceil(colors.length/(maxval-minval));
  console.log(height+(cellSize*2));
  legend.append("rect")
    .attr("x", function(d, i) { return legendElementWidth * i; })
    .attr("y", height+(cellSize*2))
    .attr("width", legendElementWidth)
    .attr("height", cellSize)
    .style("fill", function(d, i) { return colors[i*color_factor]; });

  legend.append("text")
    .attr("class", "mono")
    .text(function(d) { return d; })
    .attr("width", legendElementWidth)
    .attr("x", function(d, i) { return legendElementWidth * i; })
    .attr("y", height + (cellSize*4));

// Change ordering of cells

  function sortbylabel(rORc,i,sortOrder){
       var t = svg.transition().duration(3000);
       var log2r=[];
       var sorted; // sorted is zero-based index
       d3.selectAll(".c"+rORc+i)
         .filter(function(ce){
            log2r.push(ce.value);
          })
       ;
       if(rORc=="r"){ // sort log2ratio of a gene
         sorted=d3.range(col_number).sort(function(a,b){ if(sortOrder){ return log2r[b]-log2r[a];}else{ return log2r[a]-log2r[b];}});
         t.selectAll(".cell")
           .attr("x", function(d) { return sorted.indexOf(d.col-1) * cellSize; })
           ;
         t.selectAll(".colLabel")
          .attr("y", function (d, i) { return sorted.indexOf(i) * cellSize; })
         ;
       }else{ // sort log2ratio of a contrast
         sorted=d3.range(row_number).sort(function(a,b){if(sortOrder){ return log2r[b]-log2r[a];}else{ return log2r[a]-log2r[b];}});
         t.selectAll(".cell")
           .attr("y", function(d) { return sorted.indexOf(d.row-1) * cellSize; })
           ;
         t.selectAll(".rowLabel")
          .attr("y", function (d, i) { return sorted.indexOf(i) * cellSize; })
         ;
       }
  }

  //d3.select("#order").on("change",function(){
    this.$('#order').on("change",function(){
    order(this.value);
  });

  function order(value){
   if(value=="hclust"){
    var t = svg.transition().duration(3000);
    t.selectAll(".cell")
      .attr("x", function(d) { return hccol.indexOf(d.col) * cellSize; })
      .attr("y", function(d) { return hcrow.indexOf(d.row) * cellSize; })
      ;

    t.selectAll(".rowLabel")
      .attr("y", function (d, i) { return hcrow.indexOf(i+1) * cellSize; })
      ;

    t.selectAll(".colLabel")
      .attr("y", function (d, i) { return hccol.indexOf(i+1) * cellSize; })
      ;

   }else if (value=="probecontrast"){
    var t = svg.transition().duration(3000);
    t.selectAll(".cell")
      .attr("x", function(d) { return (d.col - 1) * cellSize; })
      .attr("y", function(d) { return (d.row - 1) * cellSize; })
      ;

    t.selectAll(".rowLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;

    t.selectAll(".colLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;

   }else if (value=="probe"){
    var t = svg.transition().duration(3000);
    t.selectAll(".cell")
      .attr("y", function(d) { return (d.row - 1) * cellSize; })
      ;

    t.selectAll(".rowLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;
   }else if (value=="contrast"){
    var t = svg.transition().duration(3000);
    t.selectAll(".cell")
      .attr("x", function(d) { return (d.col - 1) * cellSize; })
      ;
    t.selectAll(".colLabel")
      .attr("y", function (d, i) { return i * cellSize; })
      ;
   }
  }
  //
  var sa=d3.select(".g3")
      .on("mousedown", function() {
          if( !d3.event.altKey) {
             d3.selectAll(".cell-selected").classed("cell-selected",false);
             d3.selectAll(".rowLabel").classed("text-selected",false);
             d3.selectAll(".colLabel").classed("text-selected",false);
          }
         var p = d3.mouse(this);
         sa.append("rect")
         .attr({
             rx      : 0,
             ry      : 0,
             class   : "selection",
             x       : p[0],
             y       : p[1],
             width   : 1,
             height  : 1
         })
      })
      .on("mousemove", function() {
         var s = sa.select("rect.selection");

         if(!s.empty()) {
             var p = d3.mouse(this),
                 d = {
                     x       : parseInt(s.attr("x"), 10),
                     y       : parseInt(s.attr("y"), 10),
                     width   : parseInt(s.attr("width"), 10),
                     height  : parseInt(s.attr("height"), 10)
                 },
                 move = {
                     x : p[0] - d.x,
                     y : p[1] - d.y
                 }
             ;

             if(move.x < 1 || (move.x*2<d.width)) {
                 d.x = p[0];
                 d.width -= move.x;
             } else {
                 d.width = move.x;
             }

             if(move.y < 1 || (move.y*2<d.height)) {
                 d.y = p[1];
                 d.height -= move.y;
             } else {
                 d.height = move.y;
             }
             s.attr(d);

                 // deselect all temporary selected state objects
             d3.selectAll('.cell-selection.cell-selected').classed("cell-selected", false);
             d3.selectAll(".text-selection.text-selected").classed("text-selected",false);

             d3.selectAll('.cell').filter(function(cell_d, i) {
                 if(
                     !d3.select(this).classed("cell-selected") &&
                         // inner circle inside selection frame
                     (this.x.baseVal.value)+cellSize >= d.x && (this.x.baseVal.value)<=d.x+d.width &&
                     (this.y.baseVal.value)+cellSize >= d.y && (this.y.baseVal.value)<=d.y+d.height
                 ) {

                     d3.select(this)
                     .classed("cell-selection", true)
                     .classed("cell-selected", true);

                     d3.select(".r"+(cell_d.row-1))
                     .classed("text-selection",true)
                     .classed("text-selected",true);

                     d3.select(".c"+(cell_d.col-1))
                     .classed("text-selection",true)
                     .classed("text-selected",true);
                 }
             });
         }
      })
      .on("mouseup", function() {
            // remove selection frame
         sa.selectAll("rect.selection").remove();

             // remove temporary selection marker class
         d3.selectAll('.cell-selection').classed("cell-selection", false);
         d3.selectAll(".text-selection").classed("text-selection",false);
      })
      .on("mouseout", function() {
         if(d3.event.relatedTarget.tagName=='html') {
                 // remove selection frame
             sa.selectAll("rect.selection").remove();
                 // remove temporary selection marker class
             d3.selectAll('.cell-selection').classed("cell-selection", false);
             d3.selectAll(".rowLabel").classed("text-selected",false);
             d3.selectAll(".colLabel").classed("text-selected",false);
         }
      });


        }
    });
    WidgetManager.register_widget_view("HeatmapView", HeatmapView);
});
"""
