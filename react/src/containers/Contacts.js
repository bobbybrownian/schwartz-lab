import React, { Component } from 'react';
import Profile from '../components/Profile'

class Contacts extends Component {
  constructor(props){
    super(props);
    this.state = {
    profiles: [],
    pastMembers: [],
    cheese: "fondue"
    }
  }

  componentDidMount(){
    fetch(`/api/v1/profiles`)
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
        profiles: responseData.profiles,
        pastMembers: responseData.pastmembers
      })
    })
    .catch(error => console.error(`Error in fetch: ${error.message}`)); 
  }

  render(){
    let profiles = this.state.profiles.map(profile => {
      return (
        <Profile
        key={profile.id}
        id={profile.id}
        name={profile.name}
        email={profile.email}
        description={profile.description}
        img_url={profile.img_url}
        />
      )
    })

    let past_members = this.state.pastMembers.map(pastMember => {
      return (
        <p>{pastMember.name} : {pastMember.email}</p>
      )
    })
	
    return (
      <div className="contacts-header">
        <div className="row">
          <h3 className="contact-title">Principal Investigator</h3>
          <div className="columns small-4 small-offset-1">
          <img className="PI-pic" src={assetHelper["schwartz_crop2.jpg"]}></img>
          </div>
          <div className="columns small-6 PI-info">
            <ul className="PI-info">
              <li><b>Martin Alexander Schwartz</b>, Robert W. Berliner Professor of Medicine (Cardiology) and Professor of Biomedical Engineering and of Cell Biology, Yale University</li>
              <hr />
			  <li><u>email</u>: <a href="mailto:martin.schwartz@yale.edu">martin.schwartz@yale.edu</a></li>
              <li><u>phone</u>: (203) 737-2930</li>
              <li><u>address</u>: Yale Cardiovascular Research Center, 300 George Street Suite 759 New Haven, CT 06511</li>
            </ul>
          </div>
        </div>
        <div className="row members-header">
          <h3 className="contact-title">Members</h3>
          {profiles}
        </div>
        <div className="row past-member-header">
          <h3 className="contact-title">Past Members</h3>
          {past_members}
        </div>
      </div>
    )
  }
}

export default Contacts;
