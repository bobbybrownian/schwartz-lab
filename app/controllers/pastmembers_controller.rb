class PastmembersController < ApplicationController

  def new
    @pastmember = Pastmember.new
    @pastmembers = Pastmember.all
  end

  def create

    # grab params for pastmember
    pastmember_params = params["pastmember"]

    name = pastmember_params["name"]
	email = pastmember_params["email"]
    @pastmember = Pastmember.new(
                          name: name,
						  email: email
						  )
    @pastmember.save
    # binding.pry
    redirect_to maintenance_index_path
  end

  def edit

  end

  def updated

  end

  def destroy
    pastmember = Pastmember.find(params[:id])
    pastmember.destroy
    redirect_to new_pastmember_path
  end

  end
