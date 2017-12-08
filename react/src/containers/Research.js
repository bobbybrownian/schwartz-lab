import React, { Component } from 'react';
import Newstuff from '../components/Newstuff'

class Research extends Component {
  constructor(props){
    super(props)
    this.state ={
      selectedLink:"research-links",
      interests: [],
	  newstuffs: [],
      interest_topic: ""
    }
    this.getInterestLinks=this.getInterestLinks.bind(this)
  }
	
  componentDidMount(){
	fetch(`/api/v1/newstuffs`)
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
        newstuffs: responseData.newstuffs
      })
    })
    .catch(error => console.error(`Error in fetch: ${error.message}`)); 
  }
	
  getInterestLinks(id, topic){
    fetch(`/api/v1/interests/${id}`)
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
        interests: responseData,
        interest_topic:topic
      })
    })
    .catch(error => console.error(`Error in fetch: ${error.message}`));
			
	
	$('html,body').animate({
		scrollTop: $(".top-interest").offset().top},
	'slow');
	
  }
  
  render(){
    let interestsLinks = this.state.interests.map(interest =>{

      return (
        <a href={interest.link}><li>{interest.description}</li></a>
      )
    })
	
	let newstuffs = this.state.newstuffs.map(newstuff => {
      return (
        <li>{newstuff.description}</li>
      )
	})
  
    return(
      <div>
		
		<div className="row news">
			<h3 className="news-header">What's new in the Schwartz Lab?</h3>
			<ul className="news-bullet">{newstuffs}</ul>
		</div>
	  
        <div className="row top-interest">
          <div className="small-6 columns">
            <h5 className="interest-content text-center">Mechanotransduction through Integrins</h5>
            <img className="hoverImages" src={assetHelper["speckleFRET_crop.gif"]}></img>
            <span className="interest-content">
              <span className="interest-content2">
                <ul className="topic-lists">
                  <li>Understanding how cells sense applied strain and substrate stiffness through integrin-mediated adhesions.</li>
                  <li className="see-more-hover" onClick={this.getInterestLinks.bind(this, "1", "Mechanotransduction Publications")}><u>recent publications</u></li>
                </ul>
              </span>
            </span>
          </div>
          <div className="small-6 columns right-interest">
            <h5 className="interest-content text-center">Endothelial cell shear stress sensing</h5>
            <img className="hoverImages" src={assetHelper["endothelial_cells_crop.png"]}></img>
            <span className="interest-content"><span className="interest-content2">
              <ul className="topic-lists">
                <li>Understanding how endothelial cells sense shear stress from blood flow and how these pathways contribute to vascular health and disease.</li>
                <li className="see-more-hover" onClick={this.getInterestLinks.bind(this, "2", "Shear Stress Publications")}><u>recent publications</u></li>
              </ul>
            </span></span>
          </div>
		  
          
		  
		  
        </div>
        <div className="row interest-links">
          <h5><b>{this.state.interest_topic}</b></h5>
          <ol>
            {interestsLinks}
          </ol>
        </div>
		
		<div className="row interest-links2">
			<ol></ol>
		</div>

      </div>
    )
  }
}

export default Research;
