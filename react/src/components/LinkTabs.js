import React from 'react';
import { Link } from 'react-router-dom';

const LinkTabs = (props) => {
  return (
    <Link to={props.link} onClick={props.handleClick} className="labLinks">
      <div className={props.className}>
        {props.name}
      </div>
    </Link>
  )
}

export default LinkTabs
