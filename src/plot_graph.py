import numpy as np
import graph_tool.all as gt
diff_small_graph_path = "diff_graph_small_2.npz"
diff_small_graph_path = "diff_graph_small_top.npz"
def load_graph(graph_path):
    g_data = np.load(graph_path)
    labels = g_data["labels"]
    #weights = g_data["weights"]
    edges = g_data["edges"]

    g = gt.Graph(directed=False)
    g.add_vertex(len(labels))
    g.add_edge_list(edges)
    return g, labels#, weights

g, labels = load_graph(diff_small_graph_path)

genes_chromosome = np.genfromtxt("gene_chr.txt", dtype='str', usecols = [0, 1])
print(genes_chromosome.shape)

genes_to_chromosome = dict()
for i in range(genes_chromosome.shape[0]):
    genes_to_chromosome[genes_chromosome[i][0]] = int(genes_chromosome[i][1])


vchr = g.new_vertex_property("string")
vname = g.new_vertex_property("string")
vchrom = []
i = 0
for v in g.vertices():
    cur_label = labels[i].decode("utf-8")
    vchr[v] = genes_to_chromosome[cur_label]
    vname[v] = cur_label + " (" + str(genes_to_chromosome[cur_label]) + ")"
    vchrom.append(int(genes_to_chromosome[cur_label]) - 1)
    i += 1

g.vertex_properties["name"] = vname 
g.vertex_properties["chr"] = vchr 

'''
eweight = g.new_edge_property("double")
ecolor = g.new_edge_property("string")
i = 0
for e in g.edges():
    eweight[e] = weights[i]
    if weights[i] < 0:
        ecolor[e] = 'red'
    else:
        ecolor[e] = 'green'
    i += 1

g.edge_properties["weight"] = eweight
g.edge_properties["color"] = ecolor
'''

#pos = gt.planar_layout(g)
#pos = gt.radial_tree_layout(g, g.vertex(0))

for i in range(3):
    pos = gt.arf_layout(g, d = 1, a = 5, max_iter=0) # good
    #pos = gt.fruchterman_reingold_layout(g, n_iter=1000)
    #pos = gt.sfdp_layout(g, C = 1)
    #pos = gt.circle_layout(g)
    #gt.graph_draw(g, pos = pos, vertex_text=g.vertex_properties["name"], vertex_font_size=20, vertex_size=10, vertex_color = 'white', vertex_fill_color = 'blue', vertex_text_position=0, output_size=(2000, 1000), output="imgs/small_graph_top_" + str(i) + ".pdf")
    #gt.graph_draw(g, pos = pos, vertex_text=g.vertex_properties["name"], vertex_font_size=20, vertex_size=10, vertex_color = 'white', vertex_fill_color = 'blue', vertex_text_position=0, output_size=(2000, 1000), output="imgs/small_graph_top_" + str(i) + ".png")

    state = gt.minimize_blockmodel_dl(g) # , deg_corr=True, B_min = 10
    state.draw(pos=pos, vertex_shape=state.get_blocks(), vertex_text=g.vertex_properties["name"], vertex_font_size=20, vertex_size=20, edge_pen_width = 2, vertex_text_position=0, output="small_graph_top/small_graph_top_blocks_mdl_" + str(i) + ".pdf", output_size=(1500, 1000), fit_view=1.1)
    state.draw(pos=pos, vertex_shape=state.get_blocks(), vertex_text=g.vertex_properties["name"], vertex_font_size=20, vertex_size=20, edge_pen_width = 2, vertex_text_position=0, output="small_graph_top/small_graph_top_blocks_mdl_" + str(i) + ".png", output_size=(1500, 1000), fit_view=1.1)
    print(i)
    #gt.draw_hierarchy(state, layout="sfdp", vertex_text=g.vertex_properties["name"], vertex_font_size=24, vertex_text_position="centered", edge_color=g.edge_properties["color"], output_size=(2000, 1000), output="small_graph_mdl.pdf", fit_view = 0.8, hide = 2)

print(vchrom)
print(np.array(vchrom))


state = gt.NestedBlockState(g, [np.array(vchrom), np.arange(0, 22)])
gt.draw_hierarchy(state, vertex_text=g.vertex_properties["name"], vertex_font_size=24, vertex_text_position="centered", output_size=(2000, 1000), output="small_graph_top/small_graph_top_mdl.pdf", fit_view = 0.8, hide = 2)
