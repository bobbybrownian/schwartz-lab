import React, { Component } from 'react';
import PdfLinks from '../components/PdfLinks'
import SimpleMap from './SimpleMap'

class Methods extends Component {
  constructor(props){
    super(props)
    this.state ={
      selectedMethod: "",
      methodPdfs:[]
    }
    this.fetchMethodsPDF = this.fetchMethodsPDF.bind(this)
  }

  fetchMethodsPDF(method_type){
    fetch(`/api/v1/labmethods/${method_type}`)
    .then(response =>{
      if(response.ok){
        return response;
      }else {
        let errorMessage = `${response.status} (${response.statusText})`,
          error = new Error(errorMessage)
        throw(error);
      }
    })
    .then((response) => response.json())
    .then((responseData) => {
      this.setState({
        methodPdfs: responseData,
        selectedMethod:method_type
      })
    })
    .catch(error => console.error(`Error in fetch: ${error.message}`));
    }


  render(){
    let pdf_list =  this.state.methodPdfs.map(pdf=>{
      return(
        <PdfLinks
        key={pdf.id}
        id={pdf.id}
        pdf_url={pdf.pdf_url}
        title={pdf.title}
        selectedMethod={this.state.selectedMethod}
        />
      )
    })
    return(
      <div>
        <div className="row">
        </div>
        <div className="row">
          <div className="small-6 columns">
            <h2 onClick={this.fetchMethodsPDF.bind(this,"biological")} className="method-headers method-headers-top">Biological Assays</h2>
            <h2 onClick={this.fetchMethodsPDF.bind(this,"programming")} className="method-headers">Programming</h2>
            <h2 onClick={this.fetchMethodsPDF.bind(this,"philosophy")} className="method-headers">Philosophy</h2>
          </div>

		      <div className="small-6 columns method-links-layer">
            <h3>{this.state.selectedMethod}</h3>
            <ol>{pdf_list}</ol>
          </div>
        </div>
      </div>
    )
  }
}

export default Methods;
