import React from 'react'

const Profile = (props) => {
  // let methodLink = "/methods/" + props.selectedMethod + "/" + props.pdf_url
  let profilePath = "/profiles/" + props.img_url
  return (
    <div className='small-12 columns profile-section'>
      <img className="profile-picture" src={profilePath}></img>
      <div className="profile-info">
        <h3>{props.name}</h3>
        <h5>email: {props.email}</h5>
      </div>
      <p>{props.description}</p>
    </div>
  )
}

export default Profile
