import React, { Component } from 'react';
import GoogleMapReact from 'google-map-react'
import SimpleMap from './SimpleMap'

class LabHistory extends Component {
  constructor(props){
    super(props)
    this.state ={
      center: {lat: 37.09, lng: 260.71},
      zoom: 4,
	  
      selectedMethod: ""
    }
  }

  render(){
    return(
      <div>
		
        <SimpleMap />
		
      	<div className="row top-interest">
            <div className="small-7 columns historycol">
              <h5><b>A brief synopsis...</b></h5>
      		<p>My lab is among those that, in the late 1980’s, showed that integrins not only mediate physical attachment of cells to extracellular matrix but also transduce signals that inform the cell about its environment. We were responsible for identifying the Na+-H+ antiporter [1, 2], calcium channels [3, 4], c Abl [5], phosphatidylinositol 5 kinase and Rho family GTPases [6-8] as mediators of integrin signaling. We were the first to report that integrin-mediated cell adhesion promotes cell survival [9]. We were also the first to show that growth factor receptor signaling is modulated by cell adhesion and have identified several mechanisms to account for the effects [6, 10-13]. We have developed several of the widely used tools in these areas, including the pull down assay for measuring Rho activity state [14] and the fluorescence based tension sensor for measuring forces across specific proteins [15]. </p>
      		<p>A key concept throughout these studies has been that cells base their major decisions on integrated information from their environment. Our early work showed that anchorage dependence of growth, i.e., that normal cells must adhere to ECM in order to grow is because signals from integrins are required to support and sustain signals from growth factors [1, 16, 17]. Conversely, growth of cancer cells becomes anchorage-independent because downstream pathways are constitutively activated. This feature correlates closely with tumorigenesis in vivo, in particular with the ability of cancer cells to invade other tissues and metastasize to distant sites. Thus, integrin signals mediate positional specificity within multicellular organisms and constitutive activation of these pathways enable cancer cells to colonize inappropriate environments.</p>
      		<p>Later work expanded these interests into mechanotransduction, how cells sense and respond to physical stimuli. We discovered that application of fluid shear stress to vascular endothelial cells triggers conformational activation of integrins, binding of these newly activated integrins to the extracellular matrix, and subsequent signaling [18, 19] . This work led to elucidation of the junctional mechanosensory complex, the first well validated mechanotransducers for fluid shear stress [20]. It also led to identification of extracellular matrix as a key modulator of fluid shear stress signaling and atherosclerosis [21, 22]. Development of the tension sensor to measure forces across specific molecules permitted further elucidation of the molecular mechanisms in this system [23].</p>
      	  </div>
      	  <div className="small-7 columns historycol-ref">
      	  	<h5><b>References:</b></h5>
      		<ol>
      			<li>Schwartz, M.A., G. Both, and C. Lechene, The effect of cell spreading on cytoplasmic pH in normal and transformed fibroblasts. Proc. Natl. Acad. Sci. USA, 1989. 86: p. 4525-4529.</li>
      			<li>Schwartz, M.A., C. Lechene, and D.E. Ingber, Insoluble fibronectin activates the Na/H antiporter by clustering and immobilizing integrin a5b1, independent of cell shape. Proc. Natl. Acad. Sci. USA, 1991. 88: p. 7849-7853. </li>
      			<li>Schwartz, M.A., Spreading of human endothelial cells on fibronectin or vitronectin triggers elevation of intracellular free calcium. J. Cell Biol., 1993. 120: p. 1003-1010. </li>
      			<li>Schwartz, M.A., E.J. Brown, and B. Fazeli, A 50-kDa Integrin Asociated protein is required for integrin-regulated calcium entry in endothelail cells. J. Biol. Chem., 1993. 268: p. 19931-19934. </li>
      			<li>Lewis, J.M., et al., c-Abl tyrosine kinase in integrin-dependent signal transduction. Proc. Natl. Acad. Sci.USA, 1996. 93: p. 15174-15179. </li>
      			<li>McNamee, H.M., D.E. Ingber, and M.A. Schwartz, Adhesion to fibronectin stimulates inositol lipid synthesis and enhances PDGF-induced inostol lipid breakdown. J. Cell Biol., 1993. 121: p. 673-678. </li>
      			<li>Chong, L.D., et al., The small GTP-binding protein Rho regulates a phosphatidylinositol 4-phosphate 5-kinase in mammalian cells. Cell, 1994. 79: p. 507-513. </li>
      			<li>Price, L.S., et al., Activation of Rac and Cdc42 by integrins mediates cell spreading. Mol. Biol. Cell, 1998. 9: p. 1863-1871. </li>
      			<li>Meredith, J.E., B. Fazeli, and M.A. .Schwartz, The extracellular matrix as a cell survival factor. Mol. Biol. Cell, 1993. 4: p. 953-961. </li>
      			<li>Schwartz, M.A. and C. Lechene, Adhesion of 3T101/2 cells is required for protein kinase C-dependent activation of the NA-H antiporter by PDGF. Proc. Natl. Acad. Sci., 1991. 89: p. 6138-6141. </li>
      			<li>Renshaw, M.W., X.D. Ren, and M.A. Schwartz, Activation of the MAP kinase pathway by growth factors requires integrin-mediated cell adhesion. EMBO J., 1997. 16: p. 5592-5599. </li>
      			<li>delPozo, M.A., et al., Integrins regulate Rac targeting by internalization of membrane domains. Science, 2004. 303: p. 839-842. </li>
      			<li>delPozo, M.A., et al., Phospho-Caveolin-1 Mediates Integrin-Regulated Membrane Domain Internalisation. Nature Cell Biology, 2005. 7: p. 901-908. </li>
      			<li>Ren, X.D., W.B. Kiosses, and M.A. Schwartz, Regulation of the small GTP binding protein Rho by cell adhesion and the cytoskeleton. EMBO J., 1999. 18: p. 578-585. </li>
      			<li>Grashoff, C., et al., Measuring mechanical tension across vinculin reveals regulation of focal adhesion dynamics. Nature, 2010. 466(7303): p. 263-6. </li>
      			<li>Schwartz, M.A., Integrins, oncogenes and anchorage independence. J. Cell Biol., 1997. 139: p. 575-578. </li>
      			<li>Schwartz, M.A., D.Toksoz, and R. Khosravi-Far, Transformation by Rho exchange factor oncogenes is mediated by activation of an integrin-dependent pathway. EMBO J., 1996. 15: p. 6525-6530. </li>
      			<li>Tzima, E., et al., Activation of integrins in endothelial cells by fluid shear stress mediates Rho-dependent cytoskeletal alignment. EMBO J., 2001. 20: p. 4639-4647. </li>
      			<li>Tzima, E., et al., Activation of Rac in endothelial cells in response to fluid shear stress mediates gene expression and cell alignment. EMBO J., 2002. 21: p. 6791-6800. </li>
      			<li>Tzima, E., et al., Identification of a mechanosensory complex that mediates the endothelial cell response to fluid shear stress. Nature, 2005. 437: p. 426-431. </li>
      			<li>Orr, A.W., et al., The subendothelial extracellular matrix modulates NF-kappaB activation by flow: a potential role in atherosclerosis. J Cell Biol, 2005. 169(1): p. 191-202. </li>
      			<li>Orr, A.W., et al., Matrix-specific p21-activated kinase activation regulates vascular permeability in atherogenesis. J Cell Biol, 2007. 176(5): p. 719-27. </li>
      			<li>Conway, D.E., et al., Fluid shear stress on endothelial cells modulates mechanical tension across VE-cadherin and PECAM-1. Curr Biol, 2013. 23(11): p. 1024-30.</li>
      		</ol>
      	  </div>
      	  <div className="small-4 columns historycol-right">
      		<div><img src={assetHelper["speckleFRET_crop.gif"]}></img></div>
      		<div>Live imaging of NIH 3t3 cell spreading.</div>
      		<br />
      		<br />
      		<div><img src={assetHelper["orf2_crop.jpg"]}></img></div>
      		<br />
      		<div><img src={assetHelper["iki2_crop.jpg"]}></img></div>
      		<div>Integrin alpha5/alpha2 chimera knock-in mice (bottom) exhibiting less aortic atherosclerotic plaque compared to ApoE -/- mice (top). From Yun et. al. (Nature Cell Biol., 2016)</div>
      		<br />
      		<br />
      		<div><img src={assetHelper["nfkb.jpg"]}></img></div>
      		<div>Schematic of Nfkb shear pathways studied in the Schwartz Lab.</div>
      	  </div>
      	</div>
      </div>

    )
  }
}

export default LabHistory
